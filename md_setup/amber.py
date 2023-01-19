import os
import subprocess
# import tempfile
import MDAnalysis as mda

from .utils import add_hydrogen, remove_hydrogen, missing_hydrogen
from .utils import build_logger, get_formal_charge, get_n_electron
from .utils import match_pdb_to_amberff

logger = build_logger()

class AMBER_param(object):
    """
    A master function to parameterize complex pdb file with
    multiple protein chains and possibly ligands, assuming
    solvent molecules are removed from the system
    """

    def __init__(
            self, pdb,
            add_sol=True,
            keep_H=False, 
            forcefield='ff19SB', 
            forcefield_rna='OL3',
            forcefield_dna='OL15',
            watermodel='opc',
            padding=20, 
            cubic=True, 
            cation='Na+', 
            n_cations=0, 
            anion='Cl-', 
            n_anions='0'):

        self.pdb = pdb
        self.label = os.path.basename(pdb)[:-4]
        self.add_sol = add_sol
        self.keep_H = keep_H
        self.forcefield = forcefield
        self.forcefield_rna = forcefield_rna
        self.forcefield_dna = forcefield_dna
        self.watermodel = watermodel
        self.padding = padding
        self.cubic = cubic
        self.cation = cation
        self.n_cations = n_cations
        self.anion = anion 
        self.n_anions = n_anions
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
        # not_na_not_protein = mda_traj.select_atoms('not protein and not nucleic')
        # if not_na_not_protein.n_atoms == 0: 
        #     self.prot_files.append(self.pdb)
        #     return

        # loop over segments in the complex
        for seg in mda_traj.atoms.segments:
            protein = seg.atoms.select_atoms('protein')
            not_protein = seg.atoms.select_atoms('not protein')
            na_not_protein = not_protein.atoms.select_atoms('nucleic')
            not_na_not_protein = not_protein.atoms.select_atoms('not nucleic')
            print(protein.n_atoms, not_protein.n_atoms)
            # processing protein
            if protein.n_atoms != 0:
                protein_save = f'prot_seg{seg.segid}_{seg.segindex}.pdb'
                protein.write(protein_save)
                if self.keep_H and not missing_hydrogen(protein_save):
                    match_pdb_to_amberff(protein_save)
                elif not missing_hydrogen(protein_save):
                    remove_hydrogen(protein_save, protein_save)
                self.prot_files.append(protein_save)

            # processing na
            if na_not_protein.n_atoms != 0: 
                na_not_protein_save = f'na_seg{seg.segid}.pdb'
                na_not_protein.wrote(na_not_protein_save) 
                self.prot_files.append(na_not_protein_save)
            
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
        for i, lig in enumerate(self.lig_files):
            lig_tag = os.path.basename(lig)[:-4]
            if missing_hydrogen(lig):
                lig = add_hydrogen(lig)
            lig_charge = get_formal_charge(lig) if get_n_electron(lig) % 2 != 0 else 0
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
        output_top = os.path.join(output_path, f'{self.label}.prmtop')
        output_inpcrd = os.path.join(output_path, f'{self.label}.inpcrd')
        output_pdb = os.path.join(output_path, f'{self.label}.pdb')
        with open(f'leap.in', 'w+') as leap:
            # default protein ff
            leap.write(f"source leaprc.protein.{self.forcefield}\n")
            leap.write(f"source leaprc.RNA.{self.forcefield_rna}\n")
            leap.write(f"source leaprc.DNA.{self.forcefield_dna}\n")
            leap.write("source leaprc.gaff\n")
            leap.write(f"source leaprc.water.{self.watermodel}\n")
            # leap.write(f"source leaprc.")
            leap.write("set default PBRadii mbondi3\n")

            # load proteins
            prot_insts = []
            for i, prot in enumerate(self.prot_files):
                leap_name = f'biom{i}'
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
                leap.write(f"solvatebox comp {self.watermodel.upper()}BOX {self.padding}")
                if self.cubic: 
                    leap.write(f" iso\n")
                else: 
                    leap.write('/n')
                leap.write(f"addions comp {self.cation} {self.n_cations}\n")
                leap.write(f"addions comp {self.anion} {self.n_anions}\n")
            leap.write(f"saveAmberParm comp {output_top} {output_inpcrd}\n")
            leap.write(f"savepdb comp {output_pdb}\n")
            leap.write("quit\n")
        comp_info = {'top_file': output_top,
                     'inpcrd_file': output_inpcrd,
                     'pdb_file': output_pdb}
        return comp_info
