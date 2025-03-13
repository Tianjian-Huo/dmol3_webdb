import os
import re
import datetime
import sqlite3
from lib.extract_parameters import extract_parameters
from lib.save_to_db import save_to_db

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


def get_target_energy_line(energy_file):
    """解析 energy.txt，获取最后优化的 step 号及对应 init 号（若为 step 0）。"""
    with open(energy_file, 'r') as f:
        lines = f.readlines()

    last_line = lines[-1].split(":")
    last_energy = float(last_line[1].split()[1])  
    last_step_number = int(last_line[0])  

    last_matching_step = last_step_number  

    for i in range(len(lines) - 2, -1, -1):  
        parts = lines[i].split(":")
        if len(parts) > 1:
            step_energy = float(parts[1].split()[1])
            step_number = int(parts[0])

            if step_energy == last_energy:
                last_matching_step = step_number
            else:
                break  

    init_number = None
    if last_matching_step == 0:
        first_line_parts = lines[0].split(":")
        if len(first_line_parts) > 1:
            init_number = int(first_line_parts[1].split()[0])  

    return last_matching_step, init_number


def update_energy_file(energy_file, last_failed_step):
    """在 energy.txt 中回溯查找更早的 step 号，优先相同能量的 step，若无则回溯到不同能量的 step。"""
    with open(energy_file, 'r') as f:
        lines = f.readlines()

    last_energy = None
    matching_steps = []

    for line in lines:
        parts = line.split(":")
        if len(parts) > 1:
            step_number = int(parts[0])
            step_energy = float(parts[1].split()[1])

            if step_number == last_failed_step:
                last_energy = step_energy
                break

    if last_energy is None:
        return None  

    found_new_step = False
    for i in range(len(lines) - 2, -1, -1):
        parts = lines[i].split(":")
        if len(parts) > 1:
            step_number = int(parts[0])
            step_energy = float(parts[1].split()[1])

            if step_energy == last_energy and step_number < last_failed_step:
                matching_steps.append(step_number)
                found_new_step = True
            elif found_new_step:
                break  

    if matching_steps:
        return matching_steps[-1]

    for i in range(len(lines) - 2, -1, -1):
        parts = lines[i].split(":")
        if len(parts) > 1:
            step_number = int(parts[0])
            step_energy = float(parts[1].split()[1])

            if step_number < last_failed_step:
                return step_number  

    return None  


def get_step_folder(log_file, target_step, init_number=None):
    """从 log.txt 获取对应 step 的计算文件夹。"""
    with open(log_file, 'r') as f:
        lines = f.readlines()

    if target_step == 0 and init_number is not None:
        pattern = rf"init\s+{init_number}\b"
    else:
        pattern = rf"step\s+{target_step}\b"

    folder_name = None
    found_step = False

    for i, line in enumerate(lines):
        if re.search(pattern, line):
            found_step = True
            continue

        if found_step:
            if line.strip() == "" or re.match(r"step\s+\d+", line) or re.match(r"init\s+\d+", line):
                break  

            if "folder name:" in line:
                folder_name = line.split(":")[1].strip()

    return folder_name


def process_energy_file(energy_file):
    """处理 energy.txt，找到最佳 step，提取 dmol.outmol 并存入数据库。"""
    base_dir = os.path.dirname(energy_file)
    log_file = os.path.join(base_dir, "log.txt")

    if not os.path.isfile(log_file):
        return  

    last_failed_step = None

    while True:
        if last_failed_step is None:
            target_step, init_number = get_target_energy_line(energy_file)
        else:
            target_step = update_energy_file(energy_file, last_failed_step)
            init_number = None

        if target_step is None:
            return  

        if target_step == last_failed_step:
            return  

        last_failed_step = target_step  

        folder_name = get_step_folder(log_file, target_step, init_number)
        if folder_name is None:
            log_message(f"{base_dir} 未找到对应步骤的文件夹名，跳过")
            return  

        dmol_outmol_path = os.path.join(base_dir, folder_name, "dmol.outmol")

        if not os.path.isfile(dmol_outmol_path):
            continue  

        parameters, atom_species, atom_positions = extract_parameters(dmol_outmol_path)

        if parameters and atom_species and atom_positions:
            new_db_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), db_filename)
            save_to_db(new_db_path, parameters, atom_species, atom_positions)
            relative_path = os.path.relpath(base_dir, os.path.dirname(os.path.abspath(__file__)))
            log_message(f"{relative_path} 提取成功")
            return

    log_message(f"{base_dir} 提取失败，未存入数据库")


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
    root_dir = input("请输入要遍历的目录路径: ").strip()

    if not os.path.isdir(root_dir):
        print(f"错误: 目录 {root_dir} 不存在！请检查路径。")
        exit(1)

    print(f"\n开始遍历目录: {root_dir}\n")

    for dirpath, _, filenames in os.walk(root_dir):
        if "energy.txt" in filenames:
            energy_file_path = os.path.join(dirpath, "energy.txt")
            process_energy_file(energy_file_path)

    # **在 DMOL 数据提取完成前，添加一个空行**
    log_message("")  

    # **获取数据库行数**
    row_count = get_db_row_count(db_filename)
    log_message(f"DMOL 数据提取完成: {db_filename}，总行数: {row_count}")

    # **提示日志文件**
    log_message(f"日志文件已保存: {log_filename}")