import os 
import sys 

sys.path.append('..') 
from comparam.ligand import ParameterizeAMBER

pdb = os.path.abspath('./good_pl/rank113380/lig.pdb')
pro_pdb = os.path.abspath('./good_pl/rank113380/apo.pdb') 
params = ParameterizeAMBER(
    pdb, pro_pdb)

