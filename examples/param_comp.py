import os 
import sys 

sys.path.append('..') 
from comp_sim.param import ParameterizeAMBER_comp
from comp_sim.utils import to_pdb
from comp_sim.sim import simulate_explicit 

pdb = os.path.abspath('./lig1.pdb')
pro_pdb = os.path.abspath('./apo.pdb') 
params = ParameterizeAMBER_comp(
        pdb, pro_pdb)


simulate_explicit(
        params['pdb_file'], 
        params['top_file'],
        output_cm='output.h5', 
        temperature=310)



