import os
import sys
import glob
# import json
import shutil
from tqdm import tqdm
# import MDAnalysis as mda
# import pandas as pd

sys.path.append(os.path.abspath('..'))
from comp_sim.param import GMX_param

host_dir = os.getcwd()
# specify your input ligand pdbs
pdb_files = sorted(
    glob.glob('/home/hm/VM_shared/mpro/pose_san/pdbs/pdb_charmm/*po.pdb'))
print(pdb_files)
ff_dir = '/home/hm/VM_shared/mpro/pose_san/inputs/'\
        'inputs_charmm/charmm36-feb2021.ff'

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
    ff_copy = os.path.join(work_dir, os.path.basename(ff_dir))

    shutil.copy2(pdb, pdb_copy)
    print(ff_dir, ff_copy)
    shutil.copytree(ff_dir, ff_copy)
    os.chdir(work_dir)

    charmmP = GMX_param(pdb_copy, keep_H=True)
    # print(charmmP.get_box_size())
    # print(charmmP.prot_files, charmmP.lig_files)
    # charmmP.build_sys()
    os.chdir(host_dir)


    # info_list.append(info)
