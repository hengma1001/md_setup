import os
from openmm.app import PDBFile
from openmm.app.modeller import Modeller
from itertools import product
from openmm import unit

pdb_file = "ncd_dimers_0.pdb"
pdb = PDBFile(pdb_file)
label = os.path.basename(pdb_file)[:-4]
n = 2
n_copy_max = 4 # n**3
spacing = 10.0
output_filename = f"{label}_{n}_{int(spacing)}nm.pdb"

model = Modeller(pdb.topology, pdb.positions)
# model.deleteWater()
# model.delete([res for res in model.topology.residues() if "+" in res.name])
# model.delete([res for res in model.topology.residues() if "-" in res.name])
# model.delete([chain for chain in model.topology.chains() if chain.index > 1])

topology_0 = model.topology
position_lists = [model.getPositions()]
# --- For n in N ---

protein_coordinates = list(product(list(range(n)),list(range(n)),list(range(n))))[1:n_copy_max] #cubic n x n x n array of coordinates
for i, delta_xyz in enumerate(protein_coordinates):
    p_0 = position_lists[(i+1) % len(position_lists) ]
    d_p = [di* spacing for di in delta_xyz] * unit.nanometer
    position_iii = [p_0i + d_p for p_0i in p_0]
    model.add(topology_0, position_iii)

_out = open(f"{output_filename}", 'w')
PDBFile.writeFile(model.topology, model.positions, _out)
_out.close()
