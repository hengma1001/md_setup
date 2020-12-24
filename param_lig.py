import os 
import sys 
import shutil
import parmed as pmd

sys.path.append('complex_sim') 
from comp_sim.param import ParameterizeAMBER_lig

pdb = 'lig_rna_h.pdb' # 'po4.pdb'
lig_charge = -2
label = pdb.split('_')[-1][:-4]
pdb = os.path.abspath(pdb)
params = ParameterizeAMBER_lig(
    pdb, add_sol=False, lig_charge=lig_charge)
print(params) 

top_file = params['top_file']
pdb_file = params['pdb_file'] 

pdb = pmd.load_file(top_file, xyz=pdb_file) 
pdb.save(f"{label}.top") 
shutil.copy2(pdb_file, f"{label}.pdb")
