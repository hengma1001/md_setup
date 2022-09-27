import os, sys
sys.path.append("..")
from md_setup.param import NAMD_param

pdb = "pdbs/7zdq.pdb" # './spike_wt_open_prot_gl/fragments_pdb/AS1.pdb'
namd = NAMD_param(pdb, glycosylation='B') 
namd.param_comp()
