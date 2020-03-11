import os 
import subprocess
# import tempfile 
import MDAnalysis as mda 

from .utils import add_hydrogen 
from .utils import run_at_temp, to_pdb

# input_structure = './IBP.pdb' 
# pH = 7
# mol_charge = -1 


def ParameterizeAMBER_comp(pdb_lig, pdb_pro, add_sol=True):
    """
    Copied from InspireMD
    This function is pretty much a wrapper for antechamber & tleap. 
    """   
    pdb_lig = add_hydrogen(pdb_lig)
    subprocess.check_output(f'antechamber -i {pdb_lig} -fi pdb -o lig.mol2 -fo mol2 -c bcc -pf y -an y', shell=True)
    subprocess.check_output(f'parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod', shell=True)
    with open(f'leap.in', 'w+') as leap:
        leap.write("source leaprc.protein.ff14SBonlysc\n")
        leap.write("source leaprc.gaff\n")
        leap.write("source leaprc.water.tip3p\n")
        leap.write("set default PBRadii mbondi3\n")
        leap.write(f"rec = loadPDB {pdb_pro} # May need full filepath?\n")
        leap.write("saveAmberParm rec apo.prmtop apo.inpcrd\n")
        leap.write("lig = loadmol2 lig.mol2\n")
        leap.write("loadAmberParams lig.frcmod\n")
        leap.write("saveAmberParm lig lig.prmtop lig.inpcrd\n")
        leap.write("comp = combine {rec lig}\n")
        if add_sol: 
            leap.write("solvatebox comp TIP3PBOX 15\n")
            leap.write("addions comp Na+ 0\n")
            leap.write("addions comp Cl- 0\n")        
        leap.write("saveAmberParm comp comp.prmtop comp.inpcrd\n")
        leap.write("quit\n")
    subprocess.check_output(f'tleap -f leap.in', shell=True)
    if os.path.exists('comp.prmtop') and os.path.exists('comp.inpcrd'): 
        to_pdb('comp.inpcrd', 'comp.prmtop', 'comp.pdb')
        return {'top_file': os.path.abspath('comp.prmtop'), 
                'inpcrd_file': os.path.abspath('comp.inpcrd'), 
                'pdb_file': os.path.abspath('comp.pdb')}
    else: 
        raise Exception("Leap failed to build topology, check errors...")


def ParameterizeAMBER_prot(pdb_pro, add_sol=True):
    """
    This functions is to parameterize a single protein
    """
    with open('leap.in', 'w') as leap: 
        leap.write("source leaprc.protein.ff14SBonlysc\n")
        leap.write("source leaprc.gaff\n")
        leap.write("source leaprc.water.tip3p\n")
        leap.write("set default PBRadii mbondi3\n")
        leap.write(f"prot = loadPDB {pdb_pro} # May need full filepath?\n")
        if add_sol:
            leap.write("solvatebox prot TIP3PBOX 15\n")
            leap.write("addions prot Na+ 0\n")
            leap.write("addions prot Cl- 0\n")
        leap.write("saveAmberParm prot prot.prmtop prot.inpcrd\n")
        leap.write("quit\n")
    subprocess.check_output(f'tleap -f leap.in', shell=True)
    if os.path.exists('prot.prmtop') and os.path.exists('prot.inpcrd'): 
        to_pdb('prot.inpcrd', 'prot.prmtop', 'prot.pdb')
        return {'top_file': os.path.abspath('prot.prmtop'), 
                'inpcrd_file': os.path.abspath('prot.inpcrd'), 
                'pdb_file': os.path.abspath('prot.pdb')}
    else: 
        raise Exception("Leap failed to build topology, check errors...")


def ParameterizeAMBER_lig(pdb_lig, add_sol=True):
    """
    Copied from InspireMD
    This function is pretty much a wrapper for antechamber & tleap. 
    """   
    pdb_lig = add_hydrogen(pdb_lig)
    subprocess.check_output(f'antechamber -i {pdb_lig} -fi pdb -o lig.mol2 -fo mol2 -c bcc -pf y -an y', shell=True)
    subprocess.check_output(f'parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod', shell=True)
    with open(f'leap.in', 'w+') as leap:
        leap.write("source leaprc.protein.ff14SBonlysc\n")
        leap.write("source leaprc.gaff\n")
        leap.write("source leaprc.water.tip3p\n")
        leap.write("set default PBRadii mbondi3\n")
        leap.write("lig = loadmol2 lig.mol2\n")
        leap.write("loadAmberParams lig.frcmod\n")
        if add_sol: 
            leap.write("solvatebox lig TIP3PBOX 15\n")
            leap.write("addions lig Na+ 0\n")
            leap.write("addions lig Cl- 0\n")
        leap.write("saveAmberParm lig lig.prmtop lig.inpcrd\n")
        leap.write("quit\n")
    subprocess.check_output(f'tleap -f leap.in', shell=True)
    if os.path.exists('lig.prmtop') and os.path.exists('lig.inpcrd'): 
        to_pdb('lig.inpcrd', 'lig.prmtop', 'lig.pdb')
        return {'top_file': os.path.abspath('lig.prmtop'), 
                'inpcrd_file': os.path.abspath('lig.inpcrd'), 
                'pdb_file': os.path.abspath('lig.pdb')}
    else: 
        raise Exception("Leap failed to build topology, check errors...")
