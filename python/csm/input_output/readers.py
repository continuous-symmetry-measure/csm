from csm.molecule.molecule import Molecule
from csm.calculations.basic_calculations import check_perm_structure, check_perm_equivalence, check_perm_cycles
import logging
logger = logging.getLogger(__name__)




def check_perm_validity(mol, perm, **kwargs):
    falsecount, num_invalid, cycle_counts= check_perm_cycles(perm, kwargs['op_order'], kwargs['op_type'])
    if falsecount>0:
        logger.warning("Permutation does not maintain cycle structure")
    if not check_perm_equivalence(mol, perm):
        logger.warning("Permutation contains switches between non-equivalent atoms")
    try:
        if check_perm_structure(mol, perm) < 1:
            logger.warning("Permutation does not preserve molecule structure")
    except ValueError: #molecule has no structure
        pass

def read_inputs(perm_file_name=None, dir_file_name=None,  **kwargs):
    """
    Reads all the inputs for CSM
    Args:
        perm_filename: Name of permutation filename
        dir_filename: Name of direction filename
        **kwargs: Other arguments returned from get_split_arguments

    Returns:
        (molecule, perm, dir) - perm and dir may be None
    """
    molecule = Molecule.from_file(**kwargs)
    if perm_file_name:
        perm = read_perm_file(perm_file_name)
        if len(perm) != len(molecule.atoms):
            raise ValueError("Invalid permutation - permutation is of size %d but molecule has %d atoms" %
                             (len(perm), len(molecule.atoms)))
        check_perm_validity(molecule, perm, **kwargs)

    else:
        perm = None

    if dir_file_name:
        dirs = read_dir_file(dir_file_name)
        if not dirs:
            raise ValueError("provided dir file is empty")
    else:
        dirs = None
    return molecule, perm, dirs


def read_dir_file(filename):
    """
    Reads a symmetry direction file
    :param filename: Name of dir file
    :return: (x,y,z) of the symmetry axis
    """
    dirs=[]
    try:
        with open(filename, 'r') as f:
            for fileline in f:
                line = fileline.split()
                result = (float(line[0]), float(line[1]), float(line[2]))
                dirs.append(result)
        return dirs
    except:
        raise ValueError("Invalid input for use-dir")


def read_perm_file(filename):
    """
    Reads a permutation
    :param filename: Name of perm file
    :return: permutation as a list of numbers

    Check that the permutation is legal, raise ValueError if not
    """
    with open(filename, 'r') as f:
        line = f.read().split()
        used = set()

        result = []
        for num_str in line:
            try:
                num = int(num_str)
            except ValueError:
                raise ValueError("Invalid permutation %s - only numbers are allowed" % line)
            if num < 1 or num > len(line):
                raise ValueError("Invalid permutation %s - out of range" % line)
            num -= 1
            if num  in used:
                raise ValueError("Invalid permutation %s - number appears twice" % line)
            result.append(num)
            used.add(num)

    return result

