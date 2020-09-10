import os 
import re
import glob 
# import json 
import argparse 
import shutil
from tqdm import tqdm
import MDAnalysis as mda 

from comp_sim.utils import only_protein, align_to_template, clean_pdb
from comp_sim.param import ParameterizeAMBER_comp 
from comp_sim.param import ParameterizeAMBER_prot 

parser = argparse.ArgumentParser() 
parser.add_argument("-p", "--prot", 
        help="Input: protein pdb file"
        )
parser.add_argument("-l", "--lig", 
        default=None, 
        help="Input: ligand pdb file"
        ) 

args = parser.parse_args() 
host_dir = os.getcwd() 
# specify your input ligand pdbs 
pdb_file = args.lig
# specify your protein pdb 
prot_file = args.prot
# print(prot_file) 

if pdb_file == None: 
    # getting parameter for protein 
    label = 'protein' 
    pdb = prot_file
    work_dir = os.path.abspath(os.path.join(host_dir, 'input_' + label))
    os.makedirs(work_dir, exist_ok=True)
    pdb_copy = os.path.join(work_dir, os.path.basename(pdb))
    shutil.copy2(pdb, pdb_copy)
    clean_pdb(pdb_copy)
    os.chdir(work_dir) 
    info = ParameterizeAMBER_prot(pdb_copy)
    os.chdir(host_dir)
else: 
    # label for ligand identity 
    pdb = pdb_file
    pdb_code = os.path.basename(pdb).split('_')[-1][:-4] 

    work_dir = os.path.abspath(os.path.join(host_dir, 'input_' + pdb_code))
    os.makedirs(work_dir, exist_ok=True) 
    pdb_copy = os.path.join(work_dir, os.path.basename(pdb))
    prot_copy = os.path.join(work_dir, os.path.basename(prot_file)) 
    shutil.copy2(pdb, pdb_copy)
    shutil.copy2(prot_file, prot_copy) 
    clean_pdb(prot_copy) 
    clean_pdb(pdb_copy) 
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
    os.chdir(host_dir) 

