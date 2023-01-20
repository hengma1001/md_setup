import os
import h5py
import logging
import tempfile
import subprocess
import numpy as np
import MDAnalysis as mda
import parmed as pmd

# import simtk.openmm.app as app
# import simtk.openmm as omm
import simtk.unit as u

from rdkit import Chem
from mendeleev import element
from MDAnalysis.analysis import distances
from MDAnalysis.analysis import align


amber_ff_libs = [
        '/home/hm/miniconda3/envs/MD_ff/dat/leap/lib/aminoct12.lib',
        '/home/hm/miniconda3/envs/MD_ff/dat/leap/lib/aminont12.lib',
        '/home/hm/miniconda3/envs/MD_ff/dat/leap/lib/amino12.lib'
        ]

charmm_ff_libs = [
        '/home/hm/VM_shared/mpro/pose_san/inputs/'
        'inputs_charmm/charmm36-feb2021.ff/merged.rtp']
        

def build_logger(debug=0):
    logger_level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(level=logger_level, format='%(asctime)s %(message)s')
    logger = logging.getLogger(__name__)
    return logger


def trim_line(line):
    return line.split()[0].replace('"', '')


def find_diff(a, b):
    """find elements that are in A, but not in B
    """
    return sorted([i for i in a if i not in b])


def run_and_save(command, log):
    tsk = subprocess.Popen(
            command,
            stdout=log, # subprocess.PIPE,
            stderr=subprocess.STDOUT,
            shell=True)
    tsk.wait()


def get_mismatch(a, b):
    return find_diff(a, b), find_diff(b, a)


def atomList_almostEqual(a, b):
    if len(a) == len(b):
        a, b = get_atomname(a), get_atomname(b)
        # b = sorted([i[:-1] if len(i) > 1 else i for i in b])
        x = [i != j for i, j in zip(a, b)]
        if sum(x) == 0:
            return True
    return False


def get_atomname(a):
    """remove the last label in atomnames"""
    return sorted([i[:-1] if len(i) > 1 else i for i in a])


def parse_amberlib(amber_libfile):
    res_types = {}
    with open(amber_libfile, 'r') as fp:
        f_read = fp.readlines()
        line_type = None
        for line in f_read:
            if line.startswith('!!index'):
                line_type = 'index'
                continue
            elif line.startswith('!entry'):
                local_type = trim_line(line)
                if local_type.endswith('atoms'):
                    line_type = local_type.split('.')[1]
                    continue
                else:
                    line_type = None

            if line_type == 'index':
                resname = trim_line(line)
                res_types[resname] = []
            elif line_type in res_types.keys():
                res_types[line_type].append(trim_line(line))

    return res_types


def parse_charmmrtp(charmm_rtpfile):
    res_types = {}
    rtp_sections = ['bondedtypes', 'atoms', 'bonds', 'impropers', 'cmap']
    with open(charmm_rtpfile, 'r') as fp:
        f_read = fp.readlines()
        line_type = None
        resname = None
        for line in f_read:
            # parse lines
            if line.startswith(';') or line.split() == []:
                continue
            elif line.startswith('[') and line_type != 'resname':
                resname = line.split()[1]
                line_type = 'resname' if resname not in rtp_sections else None
            elif line.split()[0] == '[':
                line_type = line.split()[1]

            # parse type
            if line_type == 'resname':
                res_types[resname] = []
            elif line_type == 'atoms':
                atomname = line.split()[0]
                atomname = [atomname] if atomname != '[' else []
                res_types[resname] += atomname

    return res_types


def match_pdb_to_amberff(pdb_file, ff_libs=amber_ff_libs):
    """
    to match atom names from a pdb file to Amber forcefield
    """
    # get atom names of each residue from amber lib files
    res_types = {}
    for lib in ff_libs:
        res_types.update(parse_amberlib(lib))

    # residues with diff protonation states
    proton_res = {
            'CYS': ['CYS', 'CYM', 'CYX'],
            'HIS': ['HID', 'HIE', 'HIP']
            }

    mda_traj = mda.Universe(pdb_file)
    for res in mda_traj.residues:
        # get correct resname for termini
        resname = res.resname
        if res.resindex == 0:
            resname = 'N' + resname
        elif res.resindex == mda_traj.residues.n_residues - 1:
            resname = 'C' + resname

        # get correct resname for CYS and HIS
        atm_names = res.atoms.names
        if resname in proton_res.keys():
            for resn in proton_res[resname]:
                target_names = res_types[resn]
                mismatch = get_mismatch(atm_names, target_names)
                if atomList_almostEqual(*mismatch):
                    resname = resn
                    res.resname = resn
                    break

        target_names = res_types[resname]

        # get mismatching lists of atom, atom name and
        # force field atom name
        local_mis, target_mis = get_mismatch(atm_names, target_names)
        local_mis_atoms = [atom for name in local_mis
                           for atom in res.atoms if name == atom.name]
        ind_serach = get_atomname(target_mis)

        for atom in local_mis_atoms:
            if atom.name.startswith('HT'):
                target_name = atom.name.replace('HT', 'H')
                atom.name = target_name
            else:
                target_ind = ind_serach.index(atom.name[:-1])
                target_name = target_mis[target_ind]
                atom.name = target_name

    mda_traj.atoms.write(pdb_file)


def match_pdb_to_charmmff(pdb_file, ff_libs=charmm_ff_libs):
    """
    to match atom names from a pdb file to Charmm forcefield
    """
    # get atom names of each residue from charmm lib files
    res_types = {}
    for lib in ff_libs:
        res_types.update(parse_charmmrtp(lib))

    # residues with diff protonation states
    proton_res = {
            'CYS': ['CYS', 'CYS2', 'CYM'],
            'HIS': ['HSD', 'HSE', 'HSP'],
            'GLU': ['GLU', 'GLUP']}

    mda_traj = mda.Universe(pdb_file)
    for res in mda_traj.residues:
        # get correct resname for termini
        resname = res.resname
        # skip special residue
        if resname == 'CTL':
            continue

        # get correct resname for CYS and HIS
        atm_names = res.atoms.names
        if resname in proton_res.keys():
            for resn in proton_res[resname]:
                target_names = res_types[resn]
                mismatch = get_mismatch(atm_names, target_names)
                if atomList_almostEqual(*mismatch):
                    resname = resn
                    res.resname = resn
                    break

        target_names = res_types[resname]

        # get mismatching lists of atom, atom name and
        # force field atom name
        local_mis, target_mis = get_mismatch(atm_names, target_names)
        local_mis_atoms = [atom for name in local_mis
                           for atom in res.atoms if name == atom.name]

        for atom in local_mis_atoms:
            if atom.name.startswith('HT'):
                continue
            elif atom.name == 'OXT':
                continue
            else:
                target_ind = local_mis.index(atom.name)
                target_name = target_mis[target_ind]
                atom.name = target_name

    mda_traj.atoms.write(pdb_file)


def get_ligand(pdb_file):
    mda_trj = mda.Universe(pdb_file)
    lig = mda_trj.select_atoms('not protein')
    pdb_lig = os.path.abspath('./lig.pdb')
    lig.write(pdb_lig)
    return pdb_lig


def get_protein(pdb_file):
    mda_trj = mda.Universe(pdb_file)
    lig = mda_trj.select_atoms('protein')
    pdb_lig = os.path.abspath('./prot.pdb')
    lig.write(pdb_lig)
    return pdb_lig


def update_pdb_obabel(pdb_file, format='pdb'):
    """
    add correct conect info to pdb structure 
    obabel -ipdb lig.pdb -opdb >  lig_obabel.pdb
    """
    pdb_ob = pdb_file[:-4] + f'_ob.{format}'
    subprocess.check_output(
        f'obabel -ipdb {pdb_file} -o{format} >  {pdb_ob}',
        shell=True)
    return pdb_ob


def get_formal_charge(pdb_file, format='pdb'): 
    pdb_ob = update_pdb_obabel(pdb_file, format=format)
    if format == 'pdb':
        mol = Chem.MolFromPDBFile(pdb_ob)
    elif format == 'mol2': 
        mol = Chem.MolFromMol2File(pdb_ob)
    else: 
        raise Exception("Unknown format...")
    return Chem.GetFormalCharge(mol)


def get_lig_charge(pdb_file): 
    try:
        lig_charge = get_formal_charge(pdb_file, format='mol2')
    except: 
        lig_charge = get_formal_charge(pdb_file, format='pdb')
    n_electron = get_n_electron(pdb_file)
    if (n_electron % 2 == 0) & (lig_charge % 2 == 0): 
        return lig_charge 
    elif (n_electron % 2 != 0) & (lig_charge % 2 != 0):
        return lig_charge
    else:
        raise Exception(f"Number of electron {n_electron} and charge {lig_charge} "\
            f"are mismatch for ligand {os.path.abspath(pdb_file)}")


def get_n_electron(pdb_file): 
    mda_u = mda.Universe(pdb_file)
    n_ele = [element(atom.element).atomic_number for atom in mda_u.atoms]
    return sum(n_ele)


def is_protein(pdb_file):
    mda_trj = mda.Universe(pdb_file)
    not_prot = mda_trj.select_atoms('not protein')
    if not_prot.n_atoms == 0:
        return True
    else:
        return False


def run_at_temp(func):
    """
    Run functions at a temp dir
    """
    def wrapper(*args, **kwargs):
        current_dir = os.getcwd()
        temp_path = tempfile.TemporaryDirectory()
        os.chdir(temp_path.name)
        output = func(*args, **kwargs)
        os.chdir(current_dir)
        return output
    return wrapper


def clean_pdb(pdb_file):
    """
    Remove all entris in pdb files other than `ATOM` and HETATM`
    """
    with open(pdb_file, 'r') as pdb:
        pdb_atoms = [
            line for line in pdb
            if line.startswith('ATOM') or line.startswith('HETATM')]
    with open(pdb_file, 'w') as pdb:
        pdb.write(''.join(pdb_atoms))


def to_pdb(pos_file, top_file, pdb_file):
    top = pmd.load_file(top_file, xyz=pos_file)
    top.write_pdb(pdb_file)


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
    hydrogens = trj.select_atoms('name H*')
    missingH = True if hydrogens.n_atoms == 0 else False
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


def add_hydrogen(pdb_file):
    """
    add hydrogens to pdb structure if missing hydrogen atoms
    obabel -ipdb adp.pdb -h -opdb >  adph.pdb
    """
    if not missing_hydrogen(pdb_file):
        return pdb_file
    else:
        pdb_noH = pdb_file

    pdb_h = pdb_file[:-4] + '_h.pdb'
    subprocess.check_output(
        f'obabel -ipdb {pdb_noH} -h -opdb >  {pdb_h}',
        shell=True)
    clean_pdb(pdb_h)
    return pdb_h


def align_to_template(pdb_file, ref_file, pdb_output):
    """
    align frame to target
    """
    pdb = mda.Universe(pdb_file)
    ref = mda.Universe(ref_file)
    _ = align.alignto(pdb, ref, select='protein and name CA')
    pdb.atoms.write(pdb_output)


class ContactMapReporter(object):
    def __init__(self, file, reportInterval):
        self._file = h5py.File(file, 'w', libver='latest')
        self._file.swmr_mode = True
        self._out = self._file.create_dataset(
            'contact_maps', shape=(2, 0),
            maxshape=(None, None))
        self._reportInterval = reportInterval

    def __del__(self):
        self._file.close()

    def describeNextReport(self, simulation):
        steps = self._reportInterval - \
            simulation.currentStep % self._reportInterval
        return (steps, True, False, False, False, None)

    def report(self, simulation, state):
        ca_indices = []
        for atom in simulation.topology.atoms():
            if atom.name == 'CA':
                ca_indices.append(atom.index)
        positions = np.array(state.getPositions().value_in_unit(u.angstrom))
        # time = int(np.round(state.getTime().value_in_unit(u.picosecond)))
        positions_ca = positions[ca_indices].astype(np.float32)
        distance_matrix = distances.self_distance_array(positions_ca)
        contact_map = (distance_matrix < 8.0) * 1.0
        new_shape = (len(contact_map), self._out.shape[1] + 1)
        self._out.resize(new_shape)
        self._out[:, new_shape[1] - 1] = contact_map
        self._file.flush()
