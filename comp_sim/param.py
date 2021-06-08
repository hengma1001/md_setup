import os 
import subprocess
# import tempfile 
import MDAnalysis as mda 

from .utils import add_hydrogen, remove_hydrogen, missing_hydrogen
from .utils import get_ligand, get_protein, match_pdb_to_amberff
from .utils import run_at_temp, to_pdb, clean_pdb, is_protein


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
        """
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
        comp_info =  {'top_file': output_top, 
                      'inpcrd_file': output_inpcrd, 
                      'pdb_file': output_pdb}
        return comp_info
    





@run_at_temp
def ParameterizeAMBER_comp(pdb_lig, pdb_pro, add_sol=True, lig_charge=0):
    """
    Copied from InspireMD
    This function is pretty much a wrapper for antechamber & tleap. 
    Only suitable for single protein and single ligand 
    """   
    output_path = os.path.abspath(os.path.dirname(pdb_lig))
    output_top = os.path.join(output_path, 'comp.prmtop')
    output_inpcrd = os.path.join(output_path, 'comp.inpcrd')
    output_pdb = os.path.join(output_path, 'comp.pdb')

    # remove_hydrogen(pdb_pro, pdb_pro)
    # pdb_lig = trim_pdb(pdb_lig) 
    # if is_protein(pdb_lig): 
        # remove_hydrogen(pdb_lig, pdb_lig)
    # else: 
        # clean_pdb(pdb_lig) 
    try:
        # pdb_lig = add_hydrogen(pdb_lig) 
        # print("Using original hydrogens...") 
        subprocess.check_output(f'antechamber -i {pdb_lig} -fi pdb -o lig.mol2 -fo mol2 -c bcc -pf y -an y -nc {lig_charge}', shell=True)
    except: 
        print("incresing SCF tolerance...")
        antechamber_command = f"""antechamber -i {pdb_lig} -fi pdb -o lig.mol2 -fo mol2 -c bcc -pf y -an y -nc {lig_charge} -ek "qm_theory='AM1', grms_tol=0.0005, scfconv=1.d-9, ndiis_attempts=1000" """
        subprocess.check_output(antechamber_command, shell=True)
    subprocess.check_output(f'parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod', shell=True)
    # prepare leap input file 
    with open(f'leap.in', 'w+') as leap:
        leap.write("source leaprc.protein.ff14SB\n")
        leap.write("source leaprc.gaff\n")
        leap.write("source leaprc.water.tip3p\n")
        leap.write("set default PBRadii mbondi3\n")
        leap.write(f"rec = loadPDB {pdb_pro} # May need full filepath?\n")
        # leap.write("saveAmberParm rec apo.prmtop apo.inpcrd\n")
        if is_protein(pdb_lig): 
            leap.write(f"lig = loadPDB {pdb_lig}\n") 
        else: 
            leap.write("lig = loadmol2 lig.mol2\n")
            leap.write("loadAmberParams lig.frcmod\n")
        # leap.write("saveAmberParm lig lig.prmtop lig.inpcrd\n")
        leap.write("comp = combine {rec lig}\n")
        if add_sol: 
            leap.write("solvatebox comp TIP3PBOX 15\n")
            leap.write("addions comp Na+ 0\n")
            leap.write("addions comp Cl- 0\n")        
        leap.write(f"saveAmberParm comp {output_top} {output_inpcrd}\n")
        leap.write(f"savepdb comp {output_pdb}\n")
        leap.write("quit\n")
    subprocess.check_output(f'tleap -f leap.in', shell=True)
    if os.path.exists(output_top) and os.path.exists(output_inpcrd): 
        # to_pdb(output_inpcrd, output_top, output_pdb)
        comp_info =  {'top_file': output_top, 
                      'inpcrd_file': output_inpcrd, 
                      'pdb_file': output_pdb}
        return comp_info
    else: 
        raise Exception("Leap failed to build topology, check errors...")


@run_at_temp
def ParameterizeAMBER_comp2(pdb_comp, add_sol=True, lig_charge=0): 
    """
    This function is to build complex system with multiple chains and 
    ligands 
    """
    # protein operations 
    pdb_pro = get_protein(pdb_comp)
    # ligand operations 
    pdb_lig = get_ligand(pdb_comp)
    params = ParameterizeAMBER_comp(pdb_lig, pdb_pro, add_sol=True, lig_charge=lig_charge)
    return params


@run_at_temp
def ParameterizeAMBER_prot(pdb_pro, add_sol=True):
    """
    This functions is to parameterize a single protein
    """
    output_path = os.path.abspath(os.path.dirname(pdb_pro))
    output_top = os.path.join(output_path, 'prot.prmtop')
    output_inpcrd = os.path.join(output_path, 'prot.inpcrd')
    output_pdb = os.path.join(output_path, 'prot.pdb')
    with open('leap.in', 'w') as leap: 
        leap.write("source leaprc.protein.ff14SB\n")
        leap.write("source leaprc.gaff\n")
        leap.write("source leaprc.water.tip3p\n")
        leap.write("set default PBRadii mbondi3\n")
        leap.write(f"prot = loadPDB {pdb_pro} # May need full filepath?\n")
        if add_sol:
            leap.write("solvatebox prot TIP3PBOX 15\n")
            leap.write("addions prot Na+ 0\n")
            leap.write("addions prot Cl- 0\n")
        leap.write(f"saveAmberParm prot  {output_top} {output_inpcrd}\n")
        leap.write(f"savepdb prot {output_pdb}\n")
        leap.write("quit\n")
    subprocess.check_output(f'tleap -f leap.in', shell=True)
    if os.path.exists(output_top) and os.path.exists(output_inpcrd): 
        # to_pdb(output_inpcrd, output_top, output_pdb)
        prot_info =  {'top_file': output_top, 
                      'inpcrd_file': output_inpcrd, 
                      'pdb_file': output_pdb}
        return prot_info
    else: 
        raise Exception("Leap failed to build topology, check errors...")


@run_at_temp
def ParameterizeAMBER_lig(pdb_lig, lig_charge=0, add_sol=True):
    """
    Copied from InspireMD
    This function is pretty much a wrapper for antechamber & tleap. 
    """   
    output_path = os.path.abspath(os.path.dirname(pdb_lig))
    output_top = os.path.join(output_path, 'lig.prmtop')
    output_inpcrd = os.path.join(output_path, 'lig.inpcrd')
    output_pdb = os.path.join(output_path, 'lig.pdb')

    # pdb_lig = add_hydrogen(pdb_lig)
    subprocess.check_output(f'antechamber -i {pdb_lig} -fi pdb -o lig.mol2 -fo mol2 -c bcc -pf y -an y  -nc {lig_charge}', shell=True)
    subprocess.check_output(f'parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod', shell=True)
    with open(f'leap.in', 'w+') as leap:
        leap.write("source leaprc.protein.ff14SB\n")
        leap.write("source leaprc.gaff\n")
        leap.write("source leaprc.water.tip3p\n")
        leap.write("set default PBRadii mbondi3\n")
        leap.write("lig = loadmol2 lig.mol2\n")
        leap.write("loadAmberParams lig.frcmod\n")
        if add_sol: 
            leap.write("solvatebox lig TIP3PBOX 15\n")
            leap.write("addions lig Na+ 0\n")
            leap.write("addions lig Cl- 0\n")
        leap.write(f"saveAmberParm lig  {output_top} {output_inpcrd}\n")
        leap.write(f"savepdb lig {output_pdb}\n")
        leap.write("quit\n")
    subprocess.check_output(f'tleap -f leap.in', shell=True)

    if os.path.exists(output_top) and os.path.exists(output_inpcrd): 
        # to_pdb(output_inpcrd, output_top, output_pdb)
        lig_info =  {'top_file': output_top, 
                     'inpcrd_file': output_inpcrd, 
                     'pdb_file': output_pdb}
        return lig_info
    else: 
        raise Exception("Leap failed to build topology, check errors...")
