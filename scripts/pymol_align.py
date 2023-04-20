import os
import pymol
import glob

ref_pdb = ""
pdbs = glob.glob("./*.pdb")
pdbs = [ref_pdb] + pdbs

pymol_names = []
for pdb in pdbs:
    pdb_label = os.path.basename(pdb)
    pymol.cmd.load(pdb, pdb_label)
    pymol_names += [pdb_label]

for pymol_name in pymol_names[1:]:
    pymol.cmd.align(f'{pymol_name}', f'{pymol_names[0]}')
    pymol.cmd.save(f'{pymol_name}a.pdb', pymol_name)