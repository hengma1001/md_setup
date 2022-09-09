import os
import jinja2
import subprocess
# import tempfile
import numpy as np
import MDAnalysis as mda
from pathlib import Path
from pydantic import BaseModel

from MDAnalysis.analysis import distances

from .utils import build_logger

logger = build_logger()


class NAMD_param(object):
    """
    A master function to parameterize complex pdb file with
    multiple protein chains and possibly ligands, assuming
    solvent molecules are removed from the system
    """

    def __init__(
            self, pdb,
            add_sol=True,
            disu_cutoff:float=3,
            ):

        self.pdb = pdb
        self.pdb_label = os.path.basename(pdb)[:-4]
        self.add_sol = add_sol
        self.disu_cutoff =disu_cutoff
        self.build_segs()

    def build_segs(self):
        """
        function to separate pdb file into each component of
        the complex

        NOTE: it requires proper chain IDs for each segment
        """
        mda_traj = mda.Universe(self.pdb)
        self.segments = []
        # loop over segments in the complex
        for seg in mda_traj.atoms.segments:
            seg_save = f'{self.pdb_label}_seg{seg.segid}.pdb'
            seg.atoms.write(seg_save)
            self.segments.append(
                protein_seg(name=seg.segid, pdb=seg_save)
            )
        
        sulfur = mda_traj.select_atoms("resname CYS and name S*")
        n_sulfur = sulfur.n_atoms
        dist_disu = distances.self_distance_array(
                sulfur.atoms.positions, box=sulfur.dimensions,
            )
        pairs = [[i, j] for i in range(n_sulfur) for j in range(i+1, n_sulfur)]
        disu_pairs = np.array(pairs)[dist_disu < self.disu_cutoff]
        self.disu_strs = []
        for pair in disu_pairs: 
            disu_st = [f"{sulfur[i].segid}:{sulfur[i].resid}" for i in pair]
            self.disu_strs.append(' '.join(disu_st))

    def param_comp(self):
        """
        parameterize the complex
        """
        self.write_psfgen()
        subprocess.check_output(f'vmd -dispdev text -e psfgen.tcl', shell=True)
        
    def write_psfgen(self):
        """
        produce tleap input file
        """
        jj_env = jinja2.Environment(
            loader=jinja2.PackageLoader("comp_sim"),
            trim_blocks=True,
            lstrip_blocks=True,
            autoescape=False,
        )
        template = jj_env.get_template("psfgen.j2")
        scripts = template.render(
            {'sys_label': self.pdb_label, 
            'segs': self.segments, 
            'disu_pairs': self.disu_strs
            }
            )
        with open('psfgen.tcl', 'w') as f: 
            f.write(scripts)

    def build_glyco(self): 
        pass

        
class protein_seg(BaseModel): 
    name: str
    pdb: Path