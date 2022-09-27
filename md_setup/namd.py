import os
from typing import Optional
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

class protein_seg(BaseModel): 
    name: str
    pdb: Path

class glyco_site(BaseModel): 
    name: str
    site: str


class NAMD_param(object):
    """
    A master function to parameterize complex pdb file with
    multiple protein chains and possibly ligands, assuming
    solvent molecules are removed from the system
    """

    def __init__(
            self, pdb,
            ff_path: Path = '.',
            add_sol=True,
            disu_cutoff:float=3,
            glycosylation:Optional[str]=None,
            ):

        self.pdb = pdb
        self.ff_path = ff_path
        self.pdb_label = os.path.basename(pdb)[:-4]
        self.add_sol = add_sol
        self.disu_cutoff = disu_cutoff
        self.glycosylation = glycosylation
        self.build_segs()

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
            loader=jinja2.PackageLoader("md_setup"),
            trim_blocks=True,
            lstrip_blocks=True,
            autoescape=False,
        )
        template = jj_env.get_template("psfgen.j2")
        scripts = template.render(
            {'sys_label': self.pdb_label, 
            'ff_path': self.ff_path,
            'segs': self.segments, 
            'disu_pairs': self.disu_strs,
            'glyco_sites': self.glyco_sites,
            'add_sol': self.add_sol,
            }
            )
        with open('psfgen.tcl', 'w') as f: 
            f.write(scripts)

    def build_segs(self):
        """
        function to separate pdb file into each component of
        the complex

        NOTE: it requires proper chain IDs for each segment
        """
        mda_u = mda.Universe(self.pdb)
        self.segments = []
        # loop over segments in the complex
        for seg in mda_u.atoms.segments:
            seg_save = f'{self.pdb_label}_seg{seg.segid}.pdb'
            seg.atoms.write(seg_save)
            self.segments.append(
                protein_seg(name=seg.segid, pdb=seg_save)
            )
        self._build_disu(mda_u)
        if self.glycosylation:
            self._build_glyco(mda_u, self.glycosylation)
        else: 
            self.glyco_sites = []

    
    def _build_disu(self, mda_u):
        sulfur = mda_u.select_atoms("resname CYS and name S*")
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

    def _build_glyco(self, mda_u, glycosylation):
        """
        Find and add glyco to the N-glycosylation site on proteins if 
        required, Asn-X-Ser/Thr, X is any AA other than Pro. It's only 
        designed for the spike RBD sites for the moment as limited by 
        glyco type. 
        """ 
        if isinstance(glycosylation, str) and glycosylation in mda_u.segments.segids: 
            gly_chains = mda_u.select_atoms(f"segid {glycosylation}")
        else: 
            gly_chains = mda_u.atoms

        gly_sites = [res for res in gly_chains.residues 
                if res.resname == 'ASN' 
                and mda_u.residues.resnames[res.resindex+1] != 'PRO' 
                and mda_u.residues.resnames[res.resindex+2] in ['SER', 'THR']]
        
        glyco_sites = []
        for i, res in enumerate(gly_sites):
            site = f"{res.segid}:{res.resid}"
            glyco_sites.append(glyco_site(name=f"g{i}", site=site))
        self.glyco_sites = glyco_sites
        logger.debug(self.glyco_sites)

        
