import os
import datetime
import sqlite3
import numpy as np
from pymatgen.core import Molecule
from pymatgen.analysis.molecule_matcher import HungarianOrderMatcher, KabschMatcher
from lib.extract_parameters import extract_parameters
from lib.save_to_db import save_to_db
from lib.calculate_dos import plot_dos, read_eigenvalues
from collections import Counter
from pymatgen.core import Composition

# ========== 自动生成数据库和日志名 ==========
timestamp = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
db_filename = f"DMOL_RESULTS_{timestamp}.db"
log_filename = f"log_{timestamp}.log"

def log_message(message):
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    full_msg = f"[{ts}] {message}"
    print(full_msg)
    with open(log_filename, "a", encoding="utf-8") as f:
        f.write(full_msg + "\n")

def parse_recover_file(recover_path):
    structures = []
    with open(recover_path, 'r') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        if lines[i].startswith("pop"):
            match = lines[i].strip().split()
            pop_num = int(match[1])
            energy_line = lines[i+1].strip()
            energy = float(energy_line.split()[0])
            atom_species, atom_positions = [], []
            j = i + 2
            while j < len(lines) and lines[j].strip() and not lines[j].startswith("pop"):
                parts = lines[j].split()
                if len(parts) >= 4:
                    atom_species.append(parts[0])
                    atom_positions.append([float(x) for x in parts[1:4]])
                j += 1
            structures.append((pop_num, energy, atom_species, atom_positions))
            i = j
        else:
            i += 1
    return structures

def cluster_structures(structures, rmsd_cutoff=0.2):
    selected_popnums = []
    grouped = []

    for pop_num, energy, species, positions in structures:
        mol = Molecule(species, positions)
        added = False
        for group in grouped:
            ref_mol = Molecule(group[0][2], group[0][3])
            try:
                rmsd1 = HungarianOrderMatcher(ref_mol).fit(mol)[-1]
                rmsd2 = KabschMatcher(ref_mol).fit(mol)[-1]
                if min(rmsd1, rmsd2) < rmsd_cutoff:
                    group.append((pop_num, energy, species, positions))
                    added = True
                    break
            except Exception as e:
                continue
        if not added:
            grouped.append([(pop_num, energy, species, positions)])

    for group in grouped:
        best = min(group, key=lambda x: x[1])
        selected_popnums.append(best[0])
        log_message(f"🔹 保留结构: pop {best[0]}, energy: {best[1]} eV")

    return selected_popnums

def locate_folders_from_log(log_path, popnums):
    with open(log_path, 'r') as f:
        lines = f.readlines()

    folders = [None] * len(popnums)
    pop_index = {p: i for i, p in enumerate(popnums)}
    found = set()

    # Pass 1: 从后往前找 replace
    for i in range(len(lines) - 1, -1, -1):
        line = lines[i].strip()
        if line.startswith("replace"):
            parts = line.split()
            if parts and parts[-1].isdigit():
                pop = int(parts[-1])
                if pop in pop_index and folders[pop_index[pop]] is None:
                    # 找 folder name
                    for j in range(i - 1, -1, -1):
                        if lines[j].strip().startswith("folder name"):
                            folder = lines[j].strip().split(":")[-1].strip()
                            folders[pop_index[pop]] = folder
                            found.add(pop)
                            break
    # Pass 2: 补充剩下的 pop
    for i, pop in enumerate(popnums):
        if pop in found:
            continue
        for j in range(len(lines)):
            if lines[j].strip() == f"init {pop}":
                for k in range(j, len(lines)):
                    if lines[k].strip().startswith("folder name"):
                        folders[i] = lines[k].strip().split(":")[-1].strip()
                        break
                break
    return folders

def process_selected_folders(search_dir, folders):
    for i, folder in enumerate(folders):
        dmol_path = os.path.join(search_dir, folder, "dmol.outmol")
        if not os.path.exists(dmol_path):
            log_message(f"❌ 找不到文件: {dmol_path}")
            continue

        try:
            parameters, atom_species, atom_positions = extract_parameters(dmol_path)
        except Exception as e:
            log_message(f"❌ {dmol_path}: 提取参数失败，错误: {e}")
            continue

        if not (parameters and atom_species and atom_positions):
            log_message(f"⚠️ {dmol_path}: 信息不完整，跳过")
            continue

        formula = Composition(Counter(atom_species)).formula.replace(" ", "")
        filename = f"{formula}_{i+1}"
        parameters["filename"] = filename

        save_to_db(db_filename, parameters, atom_species, atom_positions)
        log_message(f"✅ 成功存入数据库: {dmol_path}，Filename: {filename}")

        eigenvalues, occupations = read_eigenvalues(dmol_path)
        if eigenvalues:
            if not os.path.exists("dmol_dos"):
                os.makedirs("dmol_dos")
            dos_output = os.path.join("dmol_dos", f"{filename}.png")
            plot_dos(dmol_path, dos_output, formula)
            log_message(f"📊 DOS 图已保存: {dos_output}")
        else:
            log_message(f"⚠️ 电子能级为空，跳过 DOS 绘制: {dmol_path}")

def main():
    root = input("请输入包含 search 目录的根路径: ").strip()
    if not os.path.isdir(root):
        print(f"❌ 路径不存在: {root}")
        return

    for dirpath, dirnames, _ in os.walk(root):
        if "search" in dirnames:
            search_path = os.path.join(dirpath, "search")
            recover = os.path.join(search_path, "recover.txt")
            log_txt = os.path.join(search_path, "log.txt")

            if not os.path.exists(recover) or not os.path.exists(log_txt):
                log_message(f"⚠️ 缺失 recover.txt 或 log.txt: {search_path}")
                continue

            log_message(f"📌 开始处理: {search_path}")
            structures = parse_recover_file(recover)
            popnums = cluster_structures(structures)
            folders = locate_folders_from_log(log_txt, popnums)
            process_selected_folders(search_path, folders)

    row_count = get_db_row_count(db_filename)
    log_message(f"\n✅ 所有任务完成，数据库共记录 {row_count} 条")
    log_message(f"📄 日志文件已保存: {log_filename}")

def get_db_row_count(db_path):
    if not os.path.isfile(db_path):
        return 0
    with sqlite3.connect(db_path) as conn:
        return conn.execute("SELECT COUNT(*) FROM systems").fetchone()[0]

if __name__ == "__main__":
    main()
