import os
import subprocess

# import tempfile
import MDAnalysis as mda

from .utils import (
    add_hydrogen,
    build_logger,
    get_lig_charge,
    get_lig_name,
    match_pdb_to_amberff,
    missing_hydrogen,
    remove_hydrogen,
)

logger = build_logger()


class AMBER_param(object):
    """
    A master function to parameterize complex pdb file with
    multiple protein chains and possibly ligands, assuming
    solvent molecules are removed from the system
    """

    def __init__(
        self,
        pdb,
        add_sol=True,
        keep_H=False,
        lig_charges={},
        lig_param_path="",
        forcefield="ff19SB",
        forcefield_rna=None,
        forcefield_dna=None,
        watermodel="opc",
        padding=20,
        cubic=True,
        cation="Na+",
        n_cations=0,
        anion="Cl-",
        n_anions="0",
    ):
        self.pdb = pdb
        self.label = os.path.basename(pdb)[:-4]
        self.add_sol = add_sol
        self.keep_H = keep_H
        self.lig_charges = lig_charges
        self.lig_param_path = lig_param_path
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
        self.ion_files = []
        # not_na_not_protein = mda_traj.select_atoms('not protein and not nucleic')
        # if not_na_not_protein.n_atoms == 0:
        #     self.prot_files.append(self.pdb)
        #     return

        # loop over segments in the complex
        for seg in mda_traj.atoms.segments:
            protein = seg.atoms.select_atoms("protein")
            not_protein = seg.atoms.select_atoms("not protein")
            na_not_protein = not_protein.atoms.select_atoms("nucleic")
            not_na_not_protein = not_protein.atoms.select_atoms("not nucleic")
            # print(protein.n_atoms, not_protein.n_atoms)
            # processing protein
            if protein.n_atoms != 0:
                protein_save = f"prot_seg{seg.segid}_{seg.segindex}.pdb"
                protein.write(protein_save)
                if self.keep_H and not missing_hydrogen(protein_save):
                    match_pdb_to_amberff(protein_save)
                elif not missing_hydrogen(protein_save):
                    remove_hydrogen(protein_save, protein_save)
                self.prot_files.append(protein_save)

            # processing na
            if na_not_protein.n_atoms != 0:
                na_not_protein_save = f"na_seg{seg.segid}.pdb"
                na_not_protein.wrote(na_not_protein_save)
                self.prot_files.append(na_not_protein_save)

            # processing ligands
            if not_protein.n_atoms != 0:
                for i, res in enumerate(not_protein.residues):
                    if res.atoms.n_atoms == 1:
                        ion_save = f"ion_seg{seg.segid}_{i}.pdb"
                        res.atoms.write(ion_save)
                        self.ion_files.append(ion_save)
                    else:
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
        for lig in self.lig_files:
            lig_tag = os.path.basename(lig)[:-4]
            lig_name = get_lig_name(lig)
            if missing_hydrogen(lig):
                lig = add_hydrogen(lig)
            if lig_name in self.lig_charges:
                lig_charge = self.lig_charges[lig_name]
            else:
                lig_charge = get_lig_charge(lig)

            antechamber_command = (
                f"antechamber -i {lig} -fi pdb -o {lig_tag}.mol2 "
                f"-fo mol2 -c bcc -pf y -an y -nc {lig_charge}"
            )
            subprocess.check_output(antechamber_command, shell=True)
            param_command = f"parmchk2 -i {lig_tag}.mol2 -f mol2 -o {lig_tag}.frcmod"
            subprocess.check_output(param_command, shell=True)

    def param_comp(self):
        """
        parameterize the complex
        """
        self.get_outputs()
        # skip if already run
        if os.path.exists(self.output_top) and os.path.exists(self.output_inpcrd):
            logger.info(
                f"Topology found, skipping building {os.path.basename(self.output_pdb)}..."
            )
            return

        if not self.lig_param_path:
            self.param_ligs()
        self.write_tleapIN()
        subprocess.check_output(f"tleap -f leap.in", shell=True)
        # checking whether tleap is done
        if os.path.exists(self.output_top) and os.path.exists(self.output_inpcrd):
            logger.info(f"Successfully built {os.path.basename(self.output_pdb)}...")
            return
        else:
            raise Exception("Leap failed to build topology, check errors...")

    def get_outputs(self):
        self.output_path = os.path.abspath(os.path.dirname(self.pdb))
        self.output_top = os.path.join(self.output_path, f"{self.label}.prmtop")
        self.output_inpcrd = os.path.join(self.output_path, f"{self.label}.inpcrd")
        self.output_pdb = os.path.join(self.output_path, f"{self.label}.pdb")

    def write_tleapIN(self):
        """
        produce tleap input file
        """
        with open(f"leap.in", "w+") as leap:
            # default protein ff
            if "old" in self.forcefield:
                leap.write(
                    f"source oldff/leaprc.{self.forcefield.replace('-old', '')}\n"
                )
            else:
                leap.write(f"source leaprc.protein.{self.forcefield}\n")
            if self.forcefield_dna is not None:
                leap.write(f"source leaprc.RNA.{self.forcefield_rna}\n")
            if self.forcefield_rna is not None:
                leap.write(f"source leaprc.DNA.{self.forcefield_dna}\n")
            leap.write("source leaprc.gaff\n")
            leap.write(f"source leaprc.water.{self.watermodel}\n")
            if len(self.ion_files) > 0:
                leap.write(f"loadAmberParams frcmod.ionslm_126_{self.watermodel}\n")

            # leap.write(f"source leaprc.")
            leap.write("set default PBRadii mbondi3\n")

            # load proteins
            prot_insts = []
            for i, prot in enumerate(self.prot_files):
                leap_name = f"biom{i}"
                leap.write(f"{leap_name} = loadPDB {prot}\n")
                prot_insts.append(leap_name)
            prot_insts = " ".join(prot_insts)
            # leap.write("saveAmberParm rec apo.prmtop apo.inpcrd\n")

            # load ligands
            lig_insts = []
            for lig in self.lig_files:
                lig_tag = os.path.basename(lig)[:-4]
                if self.lig_param_path:
                    lig_path_tag = f"{self.lig_param_path}/{lig_tag}"
                else:
                    lig_path_tag = lig_tag
                leap.write(f"{lig_tag} = loadmol2 {lig_path_tag}.mol2\n")
                leap.write(f"loadAmberParams {lig_path_tag}.frcmod\n")
                lig_insts.append(lig_tag)
            lig_insts = " ".join(lig_insts)

            # leap.write("saveAmberParm lig lig.prmtop lig.inpcrd\n")
            combine_insts = prot_insts + " " + lig_insts
            leap.write(f"comp = combine {{ {combine_insts} }}\n")
            if self.add_sol:
                leap.write(
                    f"solvatebox comp {self.watermodel.upper()}BOX {self.padding}"
                )
                if self.cubic:
                    leap.write(f" iso\n")
                else:
                    leap.write("/n")
                leap.write(f"addions comp {self.cation} {self.n_cations}\n")
                leap.write(f"addions comp {self.anion} {self.n_anions}\n")
            leap.write(f"saveAmberParm comp {self.output_top} {self.output_inpcrd}\n")
            leap.write(f"savepdb comp {self.output_pdb}\n")
            leap.write("quit\n")
        return 1
