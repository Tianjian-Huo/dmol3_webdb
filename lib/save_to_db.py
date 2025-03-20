from ase.db import connect
from ase import Atoms
import numpy as np

def save_to_db(db_path, parameters, atom_species, atom_positions):
    # **确保 atom_species 和 atom_positions 不为空**
    if not atom_species or not atom_positions:
        print("Error: Missing atomic data, cannot create Atoms object.")
        return

    # **构造 Atoms 对象**
    atoms = Atoms(symbols=atom_species, positions=np.array(atom_positions), pbc=[False, False, False])
    print(atoms)
    # **写入数据库**
    with connect(db_path) as db:
        db.write(
            atoms=atoms,  # **存入 Atoms 结构**
            key_value_pairs={k: v for k, v in parameters.items()}  # **存储能量等参数**
        )
        
