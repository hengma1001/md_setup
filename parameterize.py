import os 
import glob 
# import json 
import shutil
from tqdm import tqdm
import MDAnalysis as mda 
import pandas as pd 

from comp_sim.utils import is_protein, align_to_template, clean_pdb
from comp_sim.param import AMBER_param 

host_dir = os.getcwd() 
# specify your input ligand pdbs 
pdb_files = sorted(glob.glob('/home/hm/VM_shared/mpro/pose_san/pdbs/pdb_amber/*.pdb'))
print(pdb_files)

lig_charges = {'hl3': 0, 'mclue': +1}
# getting parameter for protein-ligand complexes 
info_list = []
for pdb in tqdm(pdb_files[2:]): 
    # label for ligand identity 
    pdb_code = os.path.basename(pdb)[:-4] 

    work_dir = os.path.abspath(os.path.join(host_dir, 'input_' + pdb_code))
    os.makedirs(work_dir, exist_ok=True) 
    # if glob.glob(work_dir + '/*prmtop') != []: 
    #     continue
    pdb_copy = os.path.join(work_dir, os.path.basename(pdb))

    shutil.copy2(pdb, pdb_copy)

    # clean_pdb(prot_copy) 
    # clean_pdb(pdb_copy) 
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

    # try: 
    #     print(f">>>>>>>>>{pdb_code} being parameterizaed with 0 e charge...")
    #     info = ParameterizeAMBER_comp(pdb_copy, prot_copy, lig_charge=0, add_sol=True) 
    # except: 
    #     try:
    #         print(f">>>>>>>>>{pdb_code} being parameterizaed with -1 e charge...")
    #         info = ParameterizeAMBER_comp(pdb_copy, prot_copy, lig_charge=-1, add_sol=True) 
    #     except: 
    #         print(f">>>>>>>>>{pdb_code} being parameterizaed with +1 e charge...")
    #         info = ParameterizeAMBER_comp(pdb_copy, prot_copy, lig_charge=+1, add_sol=True)


    # info_list.append(info)

