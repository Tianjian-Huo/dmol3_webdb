import os
import re
from lib.extract_parameters import extract_parameters
from lib.save_to_db import save_to_db

def get_target_energy_line(energy_file):
    """
    解析 energy.txt，获取最后优化的步骤编号和对应的初始编号（若为 step 0）。
    
    :param energy_file: energy.txt 文件路径
    :return: (最后匹配的 step 号, 初始编号 init_number 或 None)
    """
    with open(energy_file, 'r') as f:
        lines = f.readlines()

    # 解析最后一行，获取最终优化能量
    last_line = lines[-1].split(":")
    last_energy = float(last_line[1].split()[1])  # 提取能量
    last_step_number = int(last_line[0])  # 提取 step 号

    last_matching_step = last_step_number  # 记录最后匹配的 step

    # **向上查找，获取最后一个相同能量的 step**
    for i in range(len(lines) - 2, -1, -1):  
        parts = lines[i].split(":")
        if len(parts) > 1:
            step_energy = float(parts[1].split()[1])
            step_number = int(parts[0])

            if step_energy == last_energy:
                last_matching_step = step_number
            else:
                break  # 遇到不同能量时停止

    # **如果 step = 0，则查找 init number**
    init_number = None
    if last_matching_step == 0:
        first_line_parts = lines[0].split(":")
        if len(first_line_parts) > 1:
            init_number = int(first_line_parts[1].split()[0])  # 获取 init 号

    return last_matching_step, init_number


def update_energy_file(energy_file, last_failed_step):
    """
    在 energy.txt 中查找比 last_failed_step 更早的相同能量 step，
    若无匹配，则回溯到更早的不同能量 step。

    :param energy_file: energy.txt 文件路径
    :param last_failed_step: 上次失败的 step 号
    :return: 新的 step 号或 None（若无匹配）
    """
    with open(energy_file, 'r') as f:
        lines = f.readlines()

    last_energy = None
    matching_steps = []

    # **找到 last_failed_step 的能量**
    for line in lines:
        parts = line.split(":")
        if len(parts) > 1:
            step_number = int(parts[0])
            step_energy = float(parts[1].split()[1])

            if step_number == last_failed_step:
                last_energy = step_energy
                break

    if last_energy is None:
        return None  # 没有找到该 step 的能量

    # **向上查找相同的能量**
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
                break  # 遇到不同能量时停止

    # **优先返回最后匹配的 step**
    if matching_steps:
        return matching_steps[-1]

    # **如果没有相同能量的 step，则回溯到最近的不同能量 step**
    for i in range(len(lines) - 2, -1, -1):
        parts = lines[i].split(":")
        if len(parts) > 1:
            step_number = int(parts[0])
            step_energy = float(parts[1].split()[1])

            if step_number < last_failed_step:
                return step_number  # 返回更早的不同能量 step

    return None  # 无可用 step


def get_step_folder(log_file, target_step, init_number=None):
    """
    从 log.txt 获取对应 step 的计算文件夹。

    :param log_file: log.txt 文件路径
    :param target_step: 目标 step 号
    :param init_number: 目标 init 号（仅当 target_step=0 时有效）
    :return: 计算文件夹名称（若无匹配则返回 None）
    """
    with open(log_file, 'r') as f:
        lines = f.readlines()

    # **匹配 step 或 init**
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
                break  # 遇到空行或新 step 时停止

            if "folder name:" in line:
                folder_name = line.split(":")[1].strip()

    return folder_name


def process_energy_file(energy_file):
    """
    处理 energy.txt，找到最佳 step，提取 dmol.outmol 并存入数据库。

    :param energy_file: energy.txt 文件路径
    """
    base_dir = os.path.dirname(energy_file)
    log_file = os.path.join(base_dir, "log.txt")

    if not os.path.isfile(log_file):
        return  # 若 log.txt 缺失，则跳过

    last_failed_step = None

    while True:
        # **第一次查找 step**
        if last_failed_step is None:
            target_step, init_number = get_target_energy_line(energy_file)
        else:
            target_step = update_energy_file(energy_file, last_failed_step)
            init_number = None

        if target_step is None:
            return  # 无可用 step，跳过

        if target_step == last_failed_step:
            return  # 避免循环，跳过

        last_failed_step = target_step  # 记录失败的 step

        # **查找计算文件夹**
        folder_name = get_step_folder(log_file, target_step, init_number)
        if folder_name is None:
            return  # 未找到计算文件夹，跳过

        # **获取 dmol.outmol 文件路径**
        dmol_outmol_path = os.path.join(base_dir, folder_name, "dmol.outmol")

        if not os.path.isfile(dmol_outmol_path):
            continue  # 若文件缺失，则尝试下一个 step

        # **提取数据**
        parameters, atom_species, atom_positions = extract_parameters(dmol_outmol_path)

        if parameters and atom_species and atom_positions:
            new_db_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "DMOL_RESULTS.db")
            save_to_db(new_db_path, parameters, atom_species, atom_positions)
            relative_path = os.path.relpath(base_dir, os.path.dirname(os.path.abspath(__file__)))
            print(f"{relative_path} 提取成功")
            return


if __name__ == '__main__':
    # **用户输入路径**
    root_dir = input("请输入要遍历的目录路径: ").strip()

    # **检查路径是否存在**
    if not os.path.isdir(root_dir):
        print(f"错误: 目录 {root_dir} 不存在！请检查路径。")
        exit(1)

    print(f"\n开始遍历目录: {root_dir}\n")

    # **遍历目录，查找 energy.txt**
    for dirpath, _, filenames in os.walk(root_dir):
        if "energy.txt" in filenames:
            energy_file_path = os.path.join(dirpath, "energy.txt")
            process_energy_file(energy_file_path)
