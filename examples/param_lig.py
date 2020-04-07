import os 
import sys 

sys.path.append('..') 
from comp_sim.param import ParameterizeAMBER_lig

pdb = os.path.abspath('./lig1.pdb')
params = ParameterizeAMBER_lig(
    pdb, add_sol=False)
print(params) 

