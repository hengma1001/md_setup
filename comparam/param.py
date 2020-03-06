import os 
import subprocess
# import tempfile 
import MDAnalysis as mda 
from biobb_chemistry.babelm.babel_add_hydrogens import BabelAddHydrogens
from biobb_chemistry.babelm.babel_minimize import BabelMinimize
from biobb_chemistry.acpype.acpype_params_gmx import AcpypeParamsGMX

from .utils import missing_hydrogen, remove_hydrogen, run_at_temp

# input_structure = './IBP.pdb' 
# pH = 7
# mol_charge = -1 


@run_at_temp
def parameterize_ligand(
        input_structure, 
        pH=7, 
        mol_charge=0, 
        remove_H=True, 
        ):
    """
    A function that automatically parameterize a ligand 
    pdb file 

    Parameter
    ---------
    input_structure : str 
        Input structure file in pdb format, must be absolute path
    pH : float , default 7
        pH level for adding Hydrogen atoms 
    mol_charge : float, default 0
        Small molecule charge 
    readd_H : bool, default True 
        Whether to remove all the hydrogen atoms from original 
        structure 

    Return 
    ------
    ligand_dict : dict 
        Path information of created gro, itp and top files 
    """
    # input_structure = os.path.abspath(input_structure) 
    input_path = os.path.dirname(input_structure)
    ligandCode = os.path.basename(input_structure)[:-4] 

    # temp_path = tempfile.TemporaryDirectory() 
    # os.chdir(temp_path.name) 

    if missing_hydrogen(input_structure) or remove_H: 
        # 
        if not missing_hydrogen(input_structure): 
            pdb_noH = f'{ligandCode}.noH.pdb' 
            remove_hydrogen(input_structure, pdb_noH) 
        else: 
            pdb_noH = input_structure

        # add H atoms 
        output_babel_h = f'{ligandCode}.H.mol2'

        prop = {
            'ph': pH,
            'input_format': 'pdb',
            'output_format': 'mol2 '
        }

        BabelAddHydrogens(input_path=pdb_noH,
                          output_path=output_babel_h,
                          properties=prop).launch()        

    else: 
        output_babel_h = input_structure

    # Babel minimization for hydrogen 
    output_babel_min = f'{ligandCode}.H.min.pdb' 
    prop = {
        'method': 'sd',
        'criteria': '1e-10',
        'force_field': 'GAFF'
    }

    BabelMinimize(input_path=output_babel_h,
                  output_path=output_babel_min,
                  properties=prop).launch()

    # Create prop dict and inputs/outputs
    output_acpype_gro = os.path.join(input_path, ligandCode + '.gro') 
    output_acpype_itp = os.path.join(input_path, ligandCode + '.itp')
    output_acpype_top = os.path.join(input_path, ligandCode + '.top')
    output_acpype = ligandCode + 'params'
    prop = {
        'basename': output_acpype,
        'charge': mol_charge
    }

    # Create and launch bb
    AcpypeParamsGMX(input_path=output_babel_min,
                    output_path_gro=output_acpype_gro,
                    output_path_itp=output_acpype_itp,
                    output_path_top=output_acpype_top,
                    properties=prop).launch() 

    # os.chdir(input_path) 
    # temp_path.cleanup() 

    ligand_dict = {
        'name': ligandCode, 
        'pos': output_acpype_gro, 
        'itp': output_acpype_itp,  
        'top': output_acpype_top
        }
    return ligand_dict


def ParameterizeAMBER(pdb_lig, pdb_pro):
    """
    Copied from InspireMD
    This function is pretty much a wrapper for antechamber & tleap. 
    """   

    subprocess.check_output(f'antechamber -i {pdb_lig} -fi pdb -o lig.mol2 -fo mol2 -c bcc -pf y -an y', shell=True)
    subprocess.check_output(f'parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod', shell=True)
    with open(f'leap.in', 'w+') as leap:
        leap.write("source leaprc.protein.ff14SBonlysc\n")
        leap.write("source leaprc.gaff\n")
        leap.write("set default PBRadii mbondi3\n")
        leap.write(f"rec = loadPDB {pdb_pro} # May need full filepath?\n")
        leap.write("saveAmberParm rec apo.prmtop apo.inpcrd\n")
        leap.write("lig = loadmol2 lig.mol2\n")
        leap.write("loadAmberParams lig.frcmod\n")
        leap.write("com = combine {rec lig}\n")
        leap.write("saveAmberParm lig lig.prmtop lig.inpcrd\n")
        leap.write("saveAmberParm com com.prmtop com.inpcrd\n")
        leap.write("quit\n")
    subprocess.check_output(f'tleap -f leap.in', shell=True)
    if os.path.exists('com.prmtop') and os.path.exists('com.inpcrd'): 
        return {'top': os.path.abspath('com.prmtop'), 
                'pos': os.path.abspath('com.inpcrd')}
    else: 
        raise Exception("Leap failed to build topology, check errors...")


