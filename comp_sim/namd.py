import os
import jinja2
import subprocess
# import tempfile
import numpy as np
import MDAnalysis as mda
from pydantic import BaseModel
from pathlib import Path

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
            ):

        self.pdb = pdb
        self.pdb_label = pdb[:-4]
        self.add_sol = add_sol
        self.pdb_dissect()

    def pdb_dissect(self):
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

    def param_comp(self):
        """
        parameterize the complex
        """
        comp_info = self.write_tleapIN()
        subprocess.check_output(f'tleap -f leap.in', shell=True)
        if (os.path.exists(comp_info['top_file']) and
            os.path.exists(comp_info['pdb_file'])):
            return comp_info
        else:
            raise Exception("Leap failed to build topology, check errors...")

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
            'segs': self.segments}
            )
        with open('psfgen.tcl', 'w') as f: 
            f.write(scripts)


        
class protein_seg(BaseModel): 
    name: str
    pdb: Path