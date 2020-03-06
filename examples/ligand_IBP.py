import sys 

sys.path.append('..') 
from comp_sim.ligand import parameterize_ligand

pdb = './IBP.pdb'
params = parameterize_ligand(
    pdb, mol_charge=-1)

