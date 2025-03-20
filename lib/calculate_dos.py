import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import os
from scipy import constants


def read_eigenvalues(dmol_outmol_path):
    """
    解析 dmol.outmol 文件，提取电子能级（eigenvalue, eV）和占据数（occupation）。
    :param dmol_outmol_path: dmol.outmol 文件路径
    :return: eigenvalues (list), occupations (list)
    """
    eigenvalues = []
    occupations = []

    try:
        with open(dmol_outmol_path, "r") as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"❌ 文件未找到: {dmol_outmol_path}")
        return [], []

    # **反向查找电子能级部分**
    start_index = None
    for i in range(len(lines) - 1, -1, -1):
        if "state                         eigenvalue        occupation" in lines[i]:
            start_index = i + 3
            break

    if start_index is None:
        print(f"⚠️ {dmol_outmol_path}: 未找到电子能级部分")
        return [], []

    # **解析数据**
    print(f"🔍 解析电子能级: {dmol_outmol_path}")
    for line in lines[start_index:]:
        if line.strip() == "":  # 遇到空行，停止解析
            break
        parts = line.split()
        if len(parts) >= 6:
            try:
                eigenvalue_ev = float(parts[5])  # **eV**
                occupation = float(parts[6])  # **occupation**
                eigenvalues.append(-eigenvalue_ev)  # **取相反数，将负坐标移到正坐标**
                occupations.append(occupation)
            except ValueError:
                continue

    if not eigenvalues:
        print(f"⚠️ {dmol_outmol_path}: 电子能级解析失败，列表为空！")
        
    print(f"✅ 解析完成: {dmol_outmol_path}, 共解析 {len(eigenvalues)} 个能级")
    return eigenvalues, occupations


def boltzmann_weight(energies, T=300):
    """
    计算 Boltzmann 权重:
    P = exp(-E / (kB * T)) / Z
    :param energies: 能量数组
    :param T: 设定的温度（默认 300K）
    :return: 归一化权重
    """
    K_B = constants.Boltzmann / constants.e  # J/K → eV/K
    
    energies = np.array(energies)
    E_min = np.min(energies)  # **找到最小值，避免指数溢出**
    energies_shifted = energies - E_min  # **对所有能量进行平移**

    weight = np.exp(-energies_shifted / (K_B * T))

    if np.any(np.isnan(weight)) or np.sum(weight) == 0:
        print("⚠️ Boltzmann 计算失败: 权重归一化失败")
        return np.ones_like(energies)  # **返回全 1，避免 NaN 问题**

    return weight / np.sum(weight)


def calculation_dos(energy_smooth, eigenvalues, sigma, occupations=None):
    """
    计算 DOS（高斯展宽）。
    :param energy_smooth: 平滑的能量网格
    :param eigenvalues: 电子能级
    :param sigma: 高斯展宽参数
    :param occupations: 可选，占据数
    :return: 计算的 DOS
    """
    dos = np.zeros_like(energy_smooth)
    for i, occ in zip(eigenvalues, occupations or [1] * len(eigenvalues)):
        dos += norm.pdf(energy_smooth, loc=i, scale=sigma) * occ  # 考虑占据数
    return dos


def plot_dos(dmol_outmol_path, output_path, sigma=0.1, temperature=300):
    """
    计算并绘制 DOS 图。
    """
    # **读取电子能级**
    eigenvalues, occupations = read_eigenvalues(dmol_outmol_path)
    
    if not eigenvalues:
        print(f"⚠️ {output_path}: 电子能级为空，跳过绘制")
        return

    print(f"📊 开始绘制 DOS: {output_path}")
    print(f"🔍 Eigenvalues Range: min={min(eigenvalues)}, max={max(eigenvalues)}")

    # **创建保存目录（如果不存在）**
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)  # ✅ 确保目录存在
        print(f"📂 创建目录: {output_dir}")

    # **限定能量范围**
    E_min, E_max = min(eigenvalues), max(eigenvalues)
    E_range = (max(E_min - 5, -20), min(E_max + 5, 10))  # **避免极端值**
    
    # **创建平滑的能量网格**
    energy_smooth = np.linspace(E_range[0], E_range[1], 2000)

    # **计算 DOS**
    dos = calculation_dos(energy_smooth, eigenvalues, sigma, occupations)
    if np.max(dos) > 1e-10:  # **防止数值太小**
        dos /= np.max(dos)
    else:
        print(f"⚠️ DOS 归一化失败: 最大值 {np.max(dos)} 太小")
        return

    # **计算 Boltzmann 权重**
    weight = boltzmann_weight(eigenvalues, temperature)

    # **加权 DOS**
    dos_weighted = np.zeros_like(dos)
    for i in range(len(dos)):
        dos_weighted[i] = dos[i] * weight[i % len(weight)]

    # **归一化**
    if np.max(dos_weighted) != 0:
        dos_weighted /= np.max(dos_weighted)
    else:
        print(f"⚠️ {output_path}: 加权 DOS 归一化失败，最大值为 0")
        return

    # **绘制 DOS**
    plt.figure(figsize=(8, 6))
    plt.plot(energy_smooth, dos_weighted, color="black", label="DOS")
    plt.xlabel("Energy (eV)")
    plt.ylabel("DOS")
    plt.title("Density of States (Weighted)")
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)  # 增加网格线
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"✅ DOS 图已保存: {output_path}")

