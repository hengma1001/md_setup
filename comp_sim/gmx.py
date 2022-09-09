import os
import subprocess
# import tempfile
import numpy as np
import MDAnalysis as mda

from .utils import remove_hydrogen, missing_hydrogen
from .utils import build_logger
from .utils import match_pdb_to_charmmff
from .utils import run_and_save

logger = build_logger()


class GMX_param(object):
    """
    A master function to parameterize complex pdb file with
    multiple protein chains and possibly ligands in the future,
    assuming all residues in the system are available from the
    force field.

    The program might be limited in certain scenario due to
    GMX topology building process' interactive nature.
    """

    def __init__(
            self, pdb,
            add_sol=True,
            lig_charge=0,
            keep_H=False,
            box_size=0, 
            box_padding=1.,
            ion_conc=0.15):

        self.pdb = pdb
        self.add_sol = add_sol
        self.lig_charge = lig_charge
        self.keep_H = keep_H
        self.box_size = box_size
        self.box_padding = box_padding
        log_file = os.path.dirname(pdb) + '/gmx.log'
        self.log = open(log_file, 'wb')
        logger.info(f"Processing {pdb}.")
        self.build_sys()
        self.log.close()

    def pdb_preprocess(self):
        """
        function to preprocess pdb file into each component of
        the complex

        NOTE: it requires proper chain IDs for each segment
        """
        if self.keep_H and not missing_hydrogen(self.pdb):
            match_pdb_to_charmmff(self.pdb)
        elif not missing_hydrogen(self.pdb):
            remove_hydrogen(self.pdb, self.pdb)

    def top_build(self):
        """
        parameterize the ligand pdb with amber antechamber
        tools
        """
        self.top = self.pdb.replace('.pdb', '.top')
        command = f'echo -n "1\n1\n" | pdb2gmx -f {self.pdb} '\
                f'-o {self.pdb} -p {self.top}'
        run_and_save(command, self.log)

    def get_box_size(self):
        mda_traj = mda.Universe(self.pdb)
        pos_min = np.min(mda_traj.atoms.positions, axis=0)
        pos_max = np.max(mda_traj.atoms.positions, axis=0)
        box_size = np.abs(pos_max - pos_min)
        return max(box_size) / 10 + self.box_padding * 2

    def build_sol(self):
        """
        add solvent molecules
        """
        # define box size
        if self.box_size: 
            box_size = self.box_size / 10
        else: 
            box_size = self.get_box_size()
        command = f'editconf -f {self.pdb} -c'\
            f' -o {self.pdb} -bt cubic -box {box_size}'
        run_and_save(command, self.log)

        # add water molecules
        command = f"genbox -cp {self.pdb}"\
            f" -cs -p {self.top} -o {self.pdb}"
        run_and_save(command, self.log)

        # add ions
        self.tpr = self.pdb.replace('.pdb', '.tpr')
        command = f"grompp -f ions.mdp -c {self.pdb}"\
            f" -p {self.top} -o {self.tpr}"
        run_and_save(command, self.log)

        command = f"echo SOL | genion -s {self.tpr}"\
            f" -o {self.pdb} -p {self.top}"\
            f" -conc 0.15 -neutral"
        run_and_save(command, self.log)

        # verification
        command = f"grompp -f ions.mdp"\
            f" -c {self.pdb} -p {self.top}p"\
            f" -o {self.tpr}"
        tsk = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            shell=True)

        os.system("rm \\#*")

    def build_sys(self):
        logger.info("preprocessing pdb file...")
        self.pdb_preprocess()
        logger.info('building topology file...')
        self.top_build()
        logger.info('adding solvent...')
        self.build_sol()


