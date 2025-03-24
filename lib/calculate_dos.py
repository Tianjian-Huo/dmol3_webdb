import os
import numpy as np
import matplotlib.pyplot as plt

def read_eigenvalues(dmol_outmol_path):
    all_eigenvalues = []
    all_occupations = []

    try:
        with open(dmol_outmol_path, "r") as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"❌ 文件未找到: {dmol_outmol_path}")
        return [], []

    start_index = None
    for i in range(len(lines) - 1, -1, -1):
        if "state                         eigenvalue        occupation" in lines[i]:
            start_index = i + 3
            break

    if start_index is None:
        print(f"⚠️ {dmol_outmol_path}: 未找到电子能级部分")
        return [], []

    for line in lines[start_index:]:
        if line.strip() == "":
            break
        parts = line.split()
        if len(parts) >= 7:
            try:
                eigenvalue_ev = float(parts[5])
                occupation = float(parts[6])
                all_eigenvalues.append(eigenvalue_ev)
                all_occupations.append(occupation)
            except ValueError:
                continue

    homo_index = None
    for i, occ in enumerate(all_occupations):
        if occ < 1.0:
            homo_index = i
            break

    if homo_index is None:
        print("⚠️ 未找到 HOMO（occupation < 1.0）")
        return [], []

    i_start = max(homo_index - 4, 0)
    i_end = min(i_start + 9, len(all_eigenvalues) - 1)  # ✅ 向下10行（包含）

    eigenvalues = all_eigenvalues[i_start:i_end + 1]
    occupations = all_occupations[i_start:i_end + 1]

    return eigenvalues, occupations



def gaussian_broadening(eigenvalues, occupations, width=0.1, resolution=0.01, energy_window=(-20, 10)):
    min_e, max_e = energy_window
    energy_grid = np.arange(min_e, max_e, resolution)
    dos = np.zeros_like(energy_grid)

    for ev, occ in zip(eigenvalues, occupations):
        if ev < min_e - 5 * width or ev > max_e + 5 * width:
            continue
        gauss = occ * np.exp(-((energy_grid - ev) ** 2) / (2 * width ** 2))
        dos += gauss

    dos /= (width * np.sqrt(2 * np.pi))  # 归一化，单位为 states/eV
    return energy_grid, dos

def plot_dos(outmol_path, save_path, formula, width=0.1, resolution=0.01):
    eigenvalues, occupations = read_eigenvalues(outmol_path)
    if not eigenvalues:
        return

    # 自动能量范围
    ev_min = min(eigenvalues)
    ev_max = max(eigenvalues)
    margin = 1.0
    energy_window = (ev_min - margin, ev_max + margin)

    energy_axis, dos = gaussian_broadening(
        eigenvalues, occupations, width, resolution, energy_window
    )

    # 绘图
    fig, ax = plt.subplots()

    # DOS 曲线
    ax.plot(energy_axis, dos, color="#1f77b4", linewidth=1.2, label="DOS")

    # 添加能级 stick lines
    stick_height = max(dos) * 0.2
    for ev, occ in zip(eigenvalues, occupations):
        ax.vlines(ev, 0, stick_height, color='salmon', linewidth=0.5, alpha=0.6)

    # 设置坐标轴标签，使用 LaTeX
    ax.set_xlabel(r"Energy (eV)", fontsize=12)
    ax.set_ylabel(r"Density of States (eV$^{-1}$)", fontsize=12)

    # 设置分子式为左上角小标题
    ax.text(0.01, 0.95, formula, transform=ax.transAxes,
            fontsize=14, fontweight='bold', va='top', ha='left')

    # 其他图形参数
    ax.tick_params(direction='in')
    ax.grid(True, linestyle=':', linewidth=0.5)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.cla()





def log_message(msg):
    print(msg)


if __name__ == "__main__":
    current_dir = os.getcwd()
    save_dir = os.path.join(current_dir, "dmol_dos")
    os.makedirs(save_dir, exist_ok=True)

    for root, dirs, files in os.walk(current_dir):
        for file in files:
            if file.endswith(".outmol"):
                file_path = os.path.join(root, file)
                eigenvalues, occupations = read_eigenvalues(file_path)
                if eigenvalues:
                    filename = os.path.splitext(file)[0]
                    dos_output = os.path.join(save_dir, f"{filename}_dos.png")
                    plot_dos(file_path, dos_output, width=0.1, energy_window=(-20, 10))
                    log_message(f"📊 DOS 图已保存: {dos_output}")
                else:
                    log_message(f"⚠️ 电子能级为空，跳过 DOS 绘制: {file_path}")
