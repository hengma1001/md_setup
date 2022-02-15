import os
import subprocess
# import tempfile
import numpy as np
import MDAnalysis as mda

from .utils import add_hydrogen, remove_hydrogen, missing_hydrogen
from .utils import get_ligand, get_protein
from .utils import match_pdb_to_charmmff, match_pdb_to_amberff
from .utils import run_and_save
from .utils import run_at_temp, is_protein


class AMBER_param(object):
    """
    A master function to parameterize complex pdb file with
    multiple protein chains and possibly ligands, assuming
    solvent molecules are removed from the system
    """

    def __init__(
            self, pdb,
            add_sol=True,
            lig_charge=0,
            keep_H=False):

        self.pdb = pdb
        self.add_sol = add_sol
        self.lig_charge = lig_charge
        self.keep_H = keep_H
        self.pdb_dissect()

    def pdb_dissect(self):
        """
        function to separate pdb file into each component of
        the complex

        NOTE: it requires proper chain IDs for each segment
        """
        mda_traj = mda.Universe(self.pdb)
        self.prot_files = []
        self.lig_files = []
        # loop over segments in the complex
        for seg in mda_traj.atoms.segments:
            protein = seg.atoms.select_atoms('protein')
            not_protein = seg.atoms.select_atoms('not protein')
            print(protein.n_atoms, not_protein.n_atoms)
            # processing protein
            if protein.n_atoms != 0:
                protein_save = f'prot_seg{seg.segid}.pdb'
                protein.write(protein_save)
                if self.keep_H and not missing_hydrogen(protein_save):
                    match_pdb_to_amberff(protein_save)
                elif not missing_hydrogen(protein_save):
                    remove_hydrogen(protein_save, protein_save)
                self.prot_files.append(protein_save)

            # processing ligands
            if not_protein.n_atoms != 0:
                for i, res in enumerate(not_protein.residues):
                    lig_save = f"lig_seg{seg.segid}_{i}.pdb"
                    res.atoms.write(lig_save)
                    self.lig_files.append(lig_save)

    def param_ligs(self):
        """
        parameterize the ligand pdb with amber antechamber
        tools

        possible improvement:
        """
        # SCF modifications
        # antechamber_command = f"""
        #     antechamber -i {pdb_lig} -fi pdb
        #     -o lig.mol2 -fo mol2 -c bcc -pf y
        #     -an y -nc {lig_charge}
        #     -ek "qm_theory='AM1', grms_tol=0.0005,
        #     scfconv=1.d-9, ndiis_attempts=1000" """
        lig_charge = self.lig_charge
        for i, lig in enumerate(self.lig_files):
            lig_tag = os.path.basename(lig)[:-4]
            if missing_hydrogen(lig):
                lig = add_hydrogen(lig)
            antechamber_command = \
                f'antechamber -i {lig} -fi pdb -o {lig_tag}.mol2 '\
                f'-fo mol2 -c bcc -pf y -an y -nc {lig_charge}'
            subprocess.check_output(antechamber_command, shell=True)
            param_command = \
                f'parmchk2 -i {lig_tag}.mol2 -f mol2 -o {lig_tag}.frcmod'
            subprocess.check_output(param_command, shell=True)

    def param_comp(self):
        """
        parameterize the complex
        """
        self.param_ligs()
        comp_info = self.write_tleapIN()
        subprocess.check_output(f'tleap -f leap.in', shell=True)
        if (os.path.exists(comp_info['top_file']) and
            os.path.exists(comp_info['pdb_file'])):
            return comp_info
        else:
            raise Exception("Leap failed to build topology, check errors...")

    def write_tleapIN(self):
        """
        produce tleap input file
        """
        output_path = os.path.abspath(os.path.dirname(self.pdb))
        output_top = os.path.join(output_path, 'comp.prmtop')
        output_inpcrd = os.path.join(output_path, 'comp.inpcrd')
        output_pdb = os.path.join(output_path, 'comp.pdb')
        with open(f'leap.in', 'w+') as leap:
            # default protein ff
            leap.write("source leaprc.protein.ff14SB\n")
            leap.write("source leaprc.gaff\n")
            leap.write("source leaprc.water.tip3p\n")
            leap.write("set default PBRadii mbondi3\n")

            # load proteins
            prot_insts = []
            for i, prot in enumerate(self.prot_files):
                leap_name = f'rec{i}'
                leap.write(f"{leap_name} = loadPDB {prot}\n")
                prot_insts.append(leap_name)
            prot_insts = ' '.join(prot_insts)
            # leap.write("saveAmberParm rec apo.prmtop apo.inpcrd\n")

            # load ligands
            lig_insts = []
            for lig in self.lig_files:
                lig_tag = os.path.basename(lig)[:-4]
                leap.write(f"{lig_tag} = loadmol2 {lig_tag}.mol2\n")
                leap.write(f"loadAmberParams {lig_tag}.frcmod\n")
                lig_insts.append(lig_tag)
            lig_insts = ' '.join(lig_insts)

            # leap.write("saveAmberParm lig lig.prmtop lig.inpcrd\n")
            combine_insts = prot_insts + ' ' + lig_insts
            leap.write(f"comp = combine {{ {combine_insts} }}\n")
            if self.add_sol:
                leap.write("solvatebox comp TIP3PBOX 15\n")
                leap.write("addions comp Na+ 0\n")
                leap.write("addions comp Cl- 0\n")
            leap.write(f"saveAmberParm comp {output_top} {output_inpcrd}\n")
            leap.write(f"savepdb comp {output_pdb}\n")
            leap.write("quit\n")
        comp_info = {'top_file': output_top,
                     'inpcrd_file': output_inpcrd,
                     'pdb_file': output_pdb}
        return comp_info


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
        command = f'echo -n "1\n1" | pdb2gmx -f {self.pdb} '\
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
        print("preprocessing pdb file...")
        self.pdb_preprocess()
        print('building topology file...')
        self.top_build()
        print('adding solvent...')
        self.build_sol()


