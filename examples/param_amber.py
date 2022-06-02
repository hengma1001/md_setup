import os
import sys
import shutil

sys.path.append('..')
from comp_sim.param import AMBER_param

host_dir = os.getcwd()
# specify your input pdbs
pdb = 'comp.pdb'
# getting parameter for protein-ligand complexes
info_list = []

# label for ligand identity
pdb_code = os.path.basename(pdb)[:-4]

# make the work dir
work_dir = os.path.abspath(os.path.join(host_dir, 'input_' + pdb_code))
os.makedirs(work_dir, exist_ok=True)

# make a copy of pdb in the new dir
pdb_copy = os.path.join(work_dir, os.path.basename(pdb))
shutil.copy2(pdb, pdb_copy)

# run and get the parameters
os.chdir(work_dir)
amberP = AMBER_param(pdb_copy)
print(amberP.prot_files, amberP.lig_files)
amberP.param_comp()
os.chdir(host_dir)


