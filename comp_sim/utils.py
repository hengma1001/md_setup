import os 
import tempfile 
import numpy as np 
import MDAnalysis as mda 


def missing_hydrogen(pdb_file): 
    """
    Check whether a pdb file contains H atoms

    Parameters
    ----------
    pdb_file : str 
        path to input pdb file 

    Returns
    -------
    missingH : bool 
        True if missing H, false otherwise 
    """
    trj = mda.Universe(pdb_file) 
    missingH = not np.any(['H' in name for name in trj.atoms.names]) 
    return missingH


def remove_hydrogen(pdb_file, pdb_noH_file): 
    """
    remove H atoms from a pdb file 

    Parameters
    ----------
    pdb_file : str 
        path to input pdb file 
    pdb_noH_file : str 
        path to write pdb file with H removed 
    """
    trj = mda.Universe(pdb_file) 
    trj_noH = trj.select_atoms('not name H* and not name h*') 
    trj_noH.write(pdb_noH_file)


def run_at_temp(func): 
    """
    Run functions at a temp dir
    """
    def wrapper(*args, **kwargs): 
        current_dir = os.getcwd()
        temp_path = tempfile.TemporaryDirectory() 
        os.chdir(temp_path.name) 
        func(*args, **kwargs)
        os.chdir(current_dir) 
    return wrapper


def to_pdb(file_pos, file_top, file_pdb): 
    trj = mda.Universe(file_top, file_pos) 
    trj.atoms.write(file_pdb)
