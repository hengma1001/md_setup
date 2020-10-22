#!/usr/bin/env python
 
# Tell PyMOL to launch quiet (-q), fullscreen (-e) and without internal GUI (-i)
import __main__
__main__.pymol_argv = [ 'pymol', '-qei' ]
  
import pymol


import glob 

lig_pdbs = glob.glob("../pdbs/lig_*.pdb") 

for lig in lig_pdbs: 

    pymol_name = os.path.basename(lig)[:-4] 
    print(pymol_name)
    pymol.cmd.load(lig, pymol_name) 
    pymol.cmd.h_add(pymol_name) 
    pymol.cmd.save(os.path.basename(lig), pymol_name) 
