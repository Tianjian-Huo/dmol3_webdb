import re
from collections import Counter

def extract_parameters(dmol_outmol_path):
    with open(dmol_outmol_path, 'r') as f:
        lines = f.readlines()

    parameters = {}

    # 反向查找最后一次优化结果
    last_opt_section_start = None
    for i in range(len(lines) - 1, -1, -1):
        if "            Total Energy           Binding E       Cnvgnce     Time   Iter" in lines[i]:
            last_opt_section_start = i
            break

    if last_opt_section_start is None:
        print(f"{dmol_outmol_path}: No geometry optimization section found.")
        return parameters

    # 提取 HOMO_DFT 和 LUMO_DFT
    for i in range(last_opt_section_start, len(lines)):
        if "Energy of Highest Occupied Molecular Orbital" in lines[i]:
            match = re.search(r"Energy of Highest Occupied Molecular Orbital:\s*(-?\d+\.\d+)", lines[i])
            if match:
                parameters['HOMO_DFT'] = float(match.group(1)) * 27.212
        if "Energy of Lowest Unoccupied Molecular Orbital" in lines[i]:
            match = re.search(r"Energy of Lowest Unoccupied Molecular Orbital:\s*(-?\d+\.\d+)", lines[i])
            if match:
                parameters['LUMO_DFT'] = float(match.group(1)) * 27.212

    if 'HOMO_DFT' in parameters and 'LUMO_DFT' in parameters:
        parameters['GAP_DFT'] = parameters['LUMO_DFT'] - parameters['HOMO_DFT']

    # 提取 TOTEN
    for i in range(last_opt_section_start, len(lines)):
        match = re.search(r"Ef\s*(-?\d+\.\d+)", lines[i])
        if match:
            parameters['TOTEN'] = float(match.group(1)) * 27.212

    # 提取 Max_Force 并转换单位 (au → eV/Å)
    for i in range(last_opt_section_start, len(lines)):
        match = re.search(r"\|\s*\|F\|max\s*\|\s*(-?\d+\.\d+E?-?\d*)", lines[i])
        if match:
            force_au = float(match.group(1))
            force_ev_ang = force_au * 51.422067  # 1 au = 51.422067 eV/Å
            parameters['Max_Force'] = force_ev_ang

    # 提取原子信息
    atom_species = []
    atom_positions = []

    for i in range(last_opt_section_start, len(lines)):
        if "Final Coordinates (Angstroms)" in lines[i]:  # 识别坐标部分
            for j in range(i + 3, len(lines)):  # 从坐标数据行开始（跳过标题部分）
                if lines[j].strip().startswith("------"):  # 终止于分割线
                    break
                atom_data = lines[j].split()
                if len(atom_data) >= 5:  # 确保行包含完整数据
                    atom_species.append(atom_data[1])  # 提取元素符号
                    x, y, z = map(float, atom_data[2:5])  # 提取坐标
                    atom_positions.append([x, y, z])
                    
    atom_count = Counter(atom_species)

    parameters['Functional'] = "PBE"

    return parameters, atom_species, atom_positions  # ⚠️ 确保返回原子信息

