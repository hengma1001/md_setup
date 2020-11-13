import os 
import re
import glob 
# import json 
import shutil
from tqdm import tqdm
import MDAnalysis as mda 
import pandas as pd 

from comp_sim.utils import is_protein, align_to_template, clean_pdb
from comp_sim.param import ParameterizeAMBER_comp 
from comp_sim.param import ParameterizeAMBER_prot 

host_dir = os.getcwd() 
# specify your input ligand pdbs 
pdb_files = sorted(glob.glob('../lig_pdb_h/lig_*.pdb'))
# print(pdb_files)

# specify your protein pdb 
prot_file = os.path.abspath('../pdbs/3CLPro_7BQY_A_1_F_01.pdb')
# print(prot_file) 

# getting parameter for protein 
# label = 'protein' 
# pdb = prot_file
# work_dir = os.path.abspath(os.path.join(host_dir, 'input_' + label))
# os.makedirs(work_dir, exist_ok=True)
# pdb_copy = os.path.join(work_dir, os.path.basename(pdb))
# shutil.copy2(pdb, pdb_copy)
# clean_pdb(pdb_copy)
# os.chdir(work_dir) 
# info = ParameterizeAMBER_prot(pdb_copy)
# os.chdir(host_dir)

# getting parameter for protein-ligand complexes 
info_list = []
for pdb in tqdm(pdb_files): 
    # label for ligand identity 
    pdb_code = os.path.basename(pdb).split('_')[-1][:-4] 

    work_dir = os.path.abspath(os.path.join(host_dir, 'input_' + pdb_code))
    os.makedirs(work_dir, exist_ok=True) 
    if glob.glob(work_dir + '/*prmtop') != []: 
        continue
    pdb_copy = os.path.join(work_dir, os.path.basename(pdb))
    prot_copy = os.path.join(work_dir, os.path.basename(prot_file)) 
    shutil.copy2(pdb, pdb_copy)
    shutil.copy2(prot_file, prot_copy) 
    # clean_pdb(prot_copy) 
    # clean_pdb(pdb_copy) 
    os.chdir(work_dir) 
    try: 
        print(f">>>>>>>>>{pdb_code} being parameterizaed with 0 e charge...")
        info = ParameterizeAMBER_comp(pdb_copy, prot_copy, lig_charge=0, add_sol=True) 
    except: 
        try:
            print(f">>>>>>>>>{pdb_code} being parameterizaed with -1 e charge...")
            info = ParameterizeAMBER_comp(pdb_copy, prot_copy, lig_charge=-1, add_sol=True) 
        except: 
            print(f">>>>>>>>>{pdb_code} being parameterizaed with +1 e charge...")
            info = ParameterizeAMBER_comp(pdb_copy, prot_copy, lig_charge=+1, add_sol=True)


    info_list.append(info)
    os.chdir(host_dir) 

# input_filepath = os.path.abspath('./input_conf.json') 
# with open(input_filepath, 'w') as input_file: 
#     json.dump(info_list, input_file)
