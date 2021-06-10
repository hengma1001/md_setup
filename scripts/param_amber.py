import os
import glob
# import json
import shutil
from tqdm import tqdm
# import MDAnalysis as mda
# import pandas as pd

from comp_sim.param import AMBER_param

host_dir = os.getcwd()
# specify your input ligand pdbs
pdb_files = sorted(
    glob.glob('/home/hm/VM_shared/mpro/pose_san/pdbs/pdb_amber/*.pdb'))
print(pdb_files)

lig_charges = {'hl3': 0, 'mclue': +1}
# getting parameter for protein-ligand complexes
info_list = []
for pdb in tqdm(pdb_files[:]):
    # label for ligand identity
    pdb_code = os.path.basename(pdb)[:-4]

    work_dir = os.path.abspath(os.path.join(host_dir, 'input_' + pdb_code))
    os.makedirs(work_dir, exist_ok=True)
    # if glob.glob(work_dir + '/*prmtop') != []:
    #     continue
    pdb_copy = os.path.join(work_dir, os.path.basename(pdb))

    shutil.copy2(pdb, pdb_copy)

    os.chdir(work_dir)
    if pdb_code in lig_charges:
        amberP = AMBER_param(pdb_copy,
                             lig_charge=lig_charges[pdb_code],
                             keep_H=True)
    else:
        amberP = AMBER_param(pdb_copy, keep_H=True)
    print(amberP.prot_files, amberP.lig_files)
    amberP.param_comp()
    os.chdir(host_dir)

    # info_list.append(info)
