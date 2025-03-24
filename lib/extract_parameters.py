import re
from collections import Counter
from pymatgen.core.structure import Molecule
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from pymatgen.core import Composition

def fix_scientific_notation(value_str):
    """ 处理科学计数法格式，确保转换为浮点数 """
    try:
        return float(value_str.replace('E', 'e'))  # 兼容 'E' → 'e'
    except ValueError:
        return None  # 解析失败返回 None

def extract_parameters(dmol_outmol_path):
    with open(dmol_outmol_path, 'r') as f:
        lines = f.readlines()

    parameters = {}

    # **反向查找最后一次优化结果**
    last_opt_section_start = None
    for i in range(len(lines) - 1, -1, -1):
        if "            Total Energy           Binding E       Cnvgnce     Time   Iter" in lines[i]:
            last_opt_section_start = i
            break

    if last_opt_section_start is None:
        print(f"❌ {dmol_outmol_path}: No geometry optimization section found.")
        return parameters

    # **提取 HOMO_DFT 和 LUMO_DFT**
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

    # **提取 Total Energy**
    for i in range(len(lines) - 1, -1, -1):  
        if "opt==" in lines[i]:  # 识别优化行
            parts = lines[i].split()
            if len(parts) >= 3:
                try:
                    parameters['TOTEN'] = float(parts[2]) * 27.212  # **Ha → eV**
                    break
                except ValueError:
                    continue

    # **提取 Max_Force 并转换单位 (au → eV/Å)**
    for i in range(last_opt_section_start, len(lines)):
        match = re.search(r"\|\s*\|F\|max\s*\|\s*(-?\d+\.\d+E?-?\d*)", lines[i])
        if match:
            force_str = fix_scientific_notation(match.group(1))  # 处理科学计数法
            if force_str is not None:
                parameters['Max_Force'] = force_str * 51.422067  # 1 au = 51.422067 eV/Å

    # **尝试提取 Final Coordinates**
    atom_species = []
    atom_positions = []

    for i in range(last_opt_section_start, len(lines)):
        if "Final Coordinates (Angstroms)" in lines[i]:  
            for j in range(i + 3, len(lines)):  
                if lines[j].strip().startswith("------"):  
                    break
                atom_data = lines[j].split()
                if len(atom_data) >= 5:
                    atom_species.append(atom_data[1])
                    x, y, z = map(float, atom_data[2:5])
                    atom_positions.append([x, y, z])

    # **如果 Final Coordinates 为空，自动回溯**
    if not atom_species:
        geo_opt_start = None
        for i in range(len(lines) - 1, -1, -1):
            if "** GEOMETRY OPTIMIZATION IN DELOCALIZED COORDINATES **" in lines[i]:
                geo_opt_start = i
                break

        if geo_opt_start:
            for i in range(geo_opt_start, len(lines)):
                if "Input Coordinates (Angstroms)" in lines[i]:  
                    for j in range(i + 3, len(lines)):  
                        if lines[j].strip().startswith("------"):  
                            break
                        atom_data = lines[j].split()
                        if len(atom_data) >= 5:
                            atom_species.append(atom_data[1])
                            x, y, z = map(float, atom_data[2:5])
                            atom_positions.append([x, y, z])

    # **如果仍然未找到坐标，报错**
    if not atom_species:
        print(f"❌ {dmol_outmol_path}: 无法找到原子坐标！")
        return parameters

    # **统计原子种类并去重**
    unique_elements = sorted(set(atom_species))  
    parameters["Composition"] = ", ".join(unique_elements)

    # **固定计算方法**
    parameters["Calculator"] = "dmol"

    mol = Molecule(atom_species, atom_positions)

    # 获取点群
    pga = PointGroupAnalyzer(mol)
    point_group = pga.sch_symbol  # ← ✅ 这个是最兼容的方式，返回类似 "D3h", "C2v" 的字符串
    # 修正某些字符错误
    if point_group == "C*v":
        point_group = "C∞v"

    parameters["Point_Group"] = point_group

    parameters['Functional'] = "PBE"

    # **统计每种元素数量**
    element_counts = Counter(atom_species)
    formula = Composition(element_counts).formula.replace(" ", "")

    # **存入参数**
    parameters["filename"] = formula

    return parameters, atom_species, atom_positions
