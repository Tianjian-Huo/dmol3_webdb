import os
import re
import datetime
import sqlite3
import numpy as np
from pymatgen.core import Molecule
from pymatgen.analysis.molecule_matcher import HungarianOrderMatcher, KabschMatcher
from lib.extract_parameters import extract_parameters
from lib.save_to_db import save_to_db
from lib.calculate_dos import plot_dos
from lib.calculate_dos import read_eigenvalues
from collections import Counter
from pymatgen.core import Composition

# 生成时间戳文件名
timestamp = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
db_filename = f"DMOL_RESULTS_{timestamp}.db"
log_filename = f"log_{timestamp}.log"

def log_message(message):
    """在控制台输出并写入日志文件，同时附加时间戳"""
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    full_message = f"[{timestamp}] {message}"
    
    print(full_message)  
    with open(log_filename, "a", encoding="utf-8") as log_file:
        log_file.write(full_message + "\n")

def get_all_search_folders(root_dir):
    """递归查找所有包含 search 目录的路径"""
    search_folders = []
    for dirpath, dirnames, _ in os.walk(root_dir):
        if "search" in dirnames:
            search_folders.append(os.path.join(dirpath, "search"))
    return search_folders

def get_all_outmol_files(search_dir):
    """
    遍历 search 目录，查找所有 dmol.outmol 文件。
    返回 (文件路径, Molecule, TOTEN) 列表
    """
    outmol_files = []
    for root, _, files in os.walk(search_dir):
        for file in files:
            if file == "dmol.outmol":
                file_path = os.path.join(root, file)
                
                # 解析 dmol.outmol 文件
                parameters, atom_species, atom_positions = extract_parameters(file_path)
                
                if parameters and "TOTEN" in parameters:
                    molecule = get_molecule_from_outmol(file_path)
                    
                    if molecule:
                        outmol_files.append((file_path, molecule, parameters["TOTEN"]))
                        log_message(f"✅ 发现 dmol.outmol: {file_path}, TOTEN: {parameters['TOTEN']} eV")

    return outmol_files

def get_molecule_from_outmol(dmol_outmol_path):
    """解析 dmol.outmol 并构建 pymatgen 的 Molecule 对象"""
    parameters, atom_species, atom_positions = extract_parameters(dmol_outmol_path)
    if not atom_species or not atom_positions:
        return None
    return Molecule(atom_species, np.array(atom_positions))

def cluster_similar_structures(outmol_files, rmsd_cutoff=0.2):
    """
    根据 RMSD 相似度筛选结构，每组相似结构中仅保留能量最低的。
    :param outmol_files: [(文件路径, Molecule, TOTEN), ...]
    :param rmsd_cutoff: 相似度阈值
    :return: 选中的结构列表 [(文件路径, Molecule, TOTEN)]
    """
    selected_structures = []
    grouped_structures = []

    for file_path, molecule, energy in outmol_files:
        added = False
        for group in grouped_structures:
            ref_molecule, ref_energy = group[0][1], group[0][2]

            # 计算相似度
            hungarian_matcher = HungarianOrderMatcher(ref_molecule)
            hungarian_rmsd = hungarian_matcher.fit(molecule)[-1]

            kabsch_matcher = KabschMatcher(ref_molecule)
            kabsch_rmsd = kabsch_matcher.fit(molecule)[-1]

            rmsd = min(hungarian_rmsd, kabsch_rmsd)

            if rmsd < rmsd_cutoff:
                group.append((file_path, molecule, energy))
                added = True
                break

        if not added:
            grouped_structures.append([(file_path, molecule, energy)])

    # 仅保留能量最低的
    for group in grouped_structures:
        best_structure = min(group, key=lambda x: x[2])
        selected_structures.append(best_structure)
        log_message(f"🔹 选择最低能量结构: {best_structure[0]}，TOTEN: {best_structure[2]} eV")

    return selected_structures

def process_search_folder(search_dir):
    """处理单个 search 目录"""
    log_message(f"📌 开始处理 {search_dir} ...")

    outmol_files = get_all_outmol_files(search_dir)

    if not outmol_files:
        log_message(f"❌ {search_dir} 未找到 dmol.outmol 文件，跳过")
        return

    selected_structures = cluster_similar_structures(outmol_files)

    for i, (file_path, molecule, energy) in enumerate(selected_structures, start=1):
        parameters, atom_species, atom_positions = extract_parameters(file_path)

        if parameters and atom_species and atom_positions:
            # **统计每种元素数量**
            element_counts = Counter(atom_species)

            # **使用 pymatgen 生成标准化化学式**
            formula = Composition(element_counts).formula.replace(" ", "")  # 去除空格

            # **生成带编号的文件名**
            filename = f"{formula}_{i}"  # 例如 'Ca2S2_1'

            # **存入参数**
            parameters["filename"] = filename

            # **保存到数据库**
            save_to_db(db_filename, parameters, atom_species, atom_positions)
            log_message(f"✅ 成功存入数据库: {file_path}，Filename: {filename}")

            # # **提取电子能级**
            # eigenvalues, occupations = read_eigenvalues(file_path)

            # # **如果电子能级数据不为空，则绘制 DOS**
            # if len(eigenvalues) > 0:
            #     dos_output_path = os.path.join("dmol_dos", f"{filename}.png")
            #     plot_dos(file_path, dos_output_path)
            #     log_message(f"📊 DOS 图已保存: {dos_output_path}")
            # else:
            #     log_message(f"❌ {file_path}: 电子能级数据为空，跳过绘制")


def get_db_row_count(db_path):
    """获取数据库中的行数"""
    if not os.path.isfile(db_path):
        return 0
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT COUNT(*) FROM systems")
    count = cursor.fetchone()[0]
    conn.close()
    return count

if __name__ == '__main__':
    root_dir = input("请输入包含 search 目录的根路径: ").strip()

    if not os.path.isdir(root_dir):
        print(f"错误: 目录 {root_dir} 不存在！请检查路径。")
        exit(1)

    print(f"\n🔍 开始遍历 {root_dir} 下的所有 search 目录...\n")

    search_folders = get_all_search_folders(root_dir)

    if not search_folders:
        print("❌ 未找到任何 search 目录，退出")
        exit(1)

    for search_dir in search_folders:
        process_search_folder(search_dir)

    # **获取数据库行数**
    row_count = get_db_row_count(db_filename)
    log_message(f"\n✅ DMOL 数据提取完成: {db_filename}，总行数: {row_count}")

    # **提示日志文件**
    log_message(f"📄 日志文件已保存: {log_filename}")
