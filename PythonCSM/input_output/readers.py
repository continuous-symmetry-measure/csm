__author__ = 'YAEL'

def read_csm_file(f):
    """
    :param f: CSM file (the file object, not the file name)
    :return: (atoms, connectivity)
    atoms - A list of atom values: (symbol, (x,y,z))
    connectivity - A list: [[atom1 connectivity] [atom2 connectivity] ...] Connectivity lists may be empty
    """
    pass


def read_dir_file(f):
    """
    Reads a symmetry direction file
    :param f: File object
    :return: (x,y,z) of the symmetry axis
    """
    pass

def read_perm_file(f):
    """
    Reads a permutation
    :param f: File object
    :return: permutation as a list of numbers

    Check that the permutation is legal, raise ValueError if not
    """
    pass

if __name__=='__main__':
    print("Testing the file reading functions")