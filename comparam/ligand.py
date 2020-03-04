import os 
import tempfile 
import MDAnalysis as mda 
from biobb_chemistry.babelm.babel_add_hydrogens import BabelAddHydrogens
from biobb_chemistry.babelm.babel_minimize import BabelMinimize
from biobb_chemistry.acpype.acpype_params_gmx import AcpypeParamsGMX

input_structure = './IBP.pdb' 
pH = 7
mol_charge = -1 


def parameterization(
        input_structure, 
        pH=7, 
        mol_charge=-1
        ):
    """
    A function that automatically parameterize a ligand 
    pdb file 

    Parameter
    ---------
    input_structure : str 
        Input structure file in pdb format 
    pH : float 
        pH level for adding Hydrogen atoms 
    mol_charge : float 
        Small molecule charge 

    Return 
    ------
    ligand_dict : dict 
        Path information of created gro, itp and top files 
    """
    input_structure = os.path.abspath(input_structure) 
    input_path = os.path.dirname(input_structure)
    ligandCode = os.path.basename(input_structure)[:-4] 

    temp_path = tempfile.TemporaryDirectory() 
    os.chdir(temp_path.name) 

    if True:  # missing_H(input_structure): 
        # add H atoms 
        output_babel_h = f'{ligandCode}.H.pdb'

        prop = {
            'ph': pH,
            'input_format': 'pdb',
            'output_format': 'pdb'
        }

        BabelAddHydrogens(input_path=input_structure,
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

    os.chdir(input_path) 
    temp_path.cleanup() 

    ligand_dict = {
        'name': ligandCode, 
        'gro': output_acpype_gro, 
        'itp': output_acpype_itp,  
        'top': output_acpype_top
        }
    return ligand_dict
