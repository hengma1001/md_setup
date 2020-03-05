import sys 

sys.path.append('..') 
from comparam.ligand import parameterization

pdb = './IBP.pdb'
params = parameterization(
    pdb, mol_charge=-1)

