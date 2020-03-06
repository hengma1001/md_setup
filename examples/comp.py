import os 
import sys 

sys.path.append('..') 
from comparam.ligand import ParameterizeAMBER

pdb = os.path.abspath('./lig.pdb')
pro_pdb = os.path.abspath('./apo.pdb') 
params = ParameterizeAMBER(
    pdb, pro_pdb)

