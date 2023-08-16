import numpy as np

from tqdm import tqdm
from mendeleev import element
from parmed.structure import Structure
from parmed.topologyobjects import Atom

from .utils import is_float

class MAEFile():
    """to parse the maestro file to parmed structure obj"""
    def __init__(self, filename) -> None:
        self.fileobj = open(filename, 'r')
        self.struct = Structure()
        self.box = []
        self.parse()

    def parse(self):
        for line in tqdm(self.fileobj):
            self._parse_line(line)

    def _parse_line(self, line):
        line_split = line.split()
        # pbc box
        if len(line_split) == 1:
            if is_float(line_split[0]):
                self.box.append(float(line_split[0]))
        # atom info
        elif len(line_split) > 10:
            if self.struct.box is None:
                assert len(self.box) == 9
                self.struct.box = [self.box[0], self.box[4], self.box[-1], 90., 90., 90.]
            self._parse_atom_line(line)
    
    def _parse_atom_line(self, line):
        line_split = line.replace('"', '').split()
        elem = [i for i in line_split[1] if not i.isdigit()][0]
        elem = element(elem)
        atom = Atom(atomic_number=elem.atomic_number, name=line_split[1],
                    mass=elem.atomic_weight, number=line_split[0])
        atom.xx, atom.xy, atom.xz = np.array(line_split[6:9], dtype=float)
        self.struct.add_atom(atom, line_split[2], line_split[5], line_split[3])
        
