import os 
import sys 

sys.path.append('..') 
from comp_sim.param import ParameterizeAMBER_prot
from comp_sim.utils import to_pdb

pdb = os.path.abspath('./lig.pdb')
pro_pdb = os.path.abspath('./apo.pdb') 
params = ParameterizeAMBER_prot(
    pro_pdb, add_sol=False)


