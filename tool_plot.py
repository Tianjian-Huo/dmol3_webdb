import sqlite3
import matplotlib.pyplot as plt
import pandas as pd
import re
from collections import Counter
from ase.db import connect
from ase import Atoms
import numpy as np




# 📌 用户输入数据库路径
db_path = input("请输入数据库文件路径: ").strip()

# 连接数据库
conn = connect(db_path)

# 📌 提取原子结构信息
formulas = []
natoms_list = []
gap_values = []

for row in conn.select():
    atoms: Atoms = row.toatoms()  # **获取原子结构信息**
    
    # **计算原子数**
    natoms = len(atoms)
    natoms_list.append(natoms)

    # **提取分子式并去除数字**
    formula_clean = "".join(sorted(set(re.sub(r"\d+", "", atoms.get_chemical_formula(mode='hill')))))
    formulas.append(formula_clean)

    # **提取 HOMO-LUMO gap**
    gap = row.get("GAP_DFT", None)
    if gap is not None:
        gap_values.append(gap)

# **统计分子式**
formula_counts = Counter(formulas)

# **获取前 18 个最常见的分子式，剩余的归入 "Others"**
most_common = formula_counts.most_common(18)
filtered_counts = dict(most_common)
others_count = sum(count for formula, count in formula_counts.items() if formula not in filtered_counts)
filtered_counts["Others"] = others_count

# **转换为 DataFrame**
df_natoms = pd.DataFrame(natoms_list, columns=["Natoms"])
df_gap = pd.DataFrame(gap_values, columns=["GAP_DFT"])

# **绘制统计图**
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# 🔹 **分子组成统计**
axes[0].bar(filtered_counts.keys(), filtered_counts.values())  
axes[0].set_xlabel("Cluster type", fontsize=22)
axes[0].set_ylabel("Counts", fontsize=22)
axes[0].set_title("Statistics of Cluster Compositions", fontsize=22)
axes[0].tick_params(axis='x', rotation=65)  # **再多旋转一些避免重叠**
axes[0].set_xlim(-1, len(filtered_counts))  # **确保 X 轴左侧不会被截断**

# **定义 bins** (0-100 每 10 一个, 110-120 为单独 bin)
bins = list(range(0, 111, 10)) + [120]

# **映射 `>110` 到 115**
natoms_list_adjusted = [n if n <= 110 else 120 for n in natoms_list]

# **绘制直方图**
axes[1].hist(natoms_list_adjusted, bins=bins, edgecolor="black", align='mid')

# **设置 x 轴刻度**
xticks = list(range(0, 111, 10)) + [">110"]
axes[1].set_xticks(list(range(0, 111, 10)) + [120])
axes[1].set_xticklabels(xticks)  # **显示 ">110"**

# **设置标题和坐标轴**
axes[1].set_xlabel("Number of atoms in a cluster", fontsize=22)
axes[1].set_ylabel("Counts", fontsize=22)
axes[1].set_title("Statistics of Cluster Size", fontsize=22)
axes[1].set_xlim(0)  # **确保 0 贴边**

# 🔹 **HOMO-LUMO gap 统计**
axes[2].hist(df_gap["GAP_DFT"], bins=30, edgecolor='black')

# **正确显示 E_g^{DFT}**
axes[2].set_xlabel(r"$\mathrm{E_g^{DFT}}$ (eV)", fontsize=22)
axes[2].set_ylabel("Counts", fontsize=22)
axes[2].set_title("Statistics of HOMO-LUMO Gap", fontsize=22)
axes[2].set_xlim(0, max(gap_values) + 0.5)  # **确保 0 贴左边**
axes[2].set_xticks(np.linspace(0, max(gap_values), num=10))  # **均匀划分 X 轴刻度**

# **调整布局**
plt.tight_layout()
plt.show()
