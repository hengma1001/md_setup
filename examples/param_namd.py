import os, sys
sys.path.append(os.path.abspath('..'))
from comp_sim.param import NAMD_param

pdb = './spike_wt_open_prot_gl/fragments_pdb/AS1.pdb'
namd = NAMD_param(pdb, glycosylation=True) 
namd.param_comp()
