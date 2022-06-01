import os
import sys
import shutil

sys.path.append(os.path.abspath('..'))
from comp_sim.param import GMX_param

host_dir = os.getcwd()
# specify your input pdb
pdb = 'apo.pdb'
ff_dir = os.path.abspath('charmm36-feb2021.ff') 

# label for ligand identity
pdb_code = os.path.basename(pdb)[:-4]

# make a work dir
work_dir = os.path.abspath(os.path.join(host_dir, 'input_' + pdb_code))
os.makedirs(work_dir, exist_ok=True)

# make a copy of pdb and ff
pdb_copy = os.path.join(work_dir, os.path.basename(pdb))
ff_copy = os.path.join(work_dir, os.path.basename(ff_dir))
shutil.copy2(pdb, pdb_copy)
shutil.copytree(ff_dir, ff_copy)

# run the parameterization
os.chdir(work_dir)
charmmP = GMX_param(pdb_copy, box_size=150)
os.chdir(host_dir)
