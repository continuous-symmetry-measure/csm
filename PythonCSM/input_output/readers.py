from openbabel import OBFormat
import openbabel

__author__ = 'YAEL'

from input_output.molecule import Atom


def open_non_csm_file(args_dict):
    """
    :param args_dict: dictionary of processed command line arguments
    :return: OBMol object created from input file by OpenBabel
    """
    conv = openbabel.OBConversion()
    mol = openbabel.OBMol()
    if not args_dict['useformat']:
        ob_format = conv.FormatFromExt(args_dict['inFileName'])
        if not ob_format:
            raise ValueError("Error discovering format from filename " + args_dict['inFileName'])
        if not conv.SetInFormat(ob_format):
            raise ValueError("Error setting openbabel format")
    else:
        if not conv.SetInFormat(args_dict['format']):
            raise ValueError("Error setting input format to " + args_dict['format'])

    if not args_dict['babelBond']:
        conv.SetOptions("b", conv.INOPTIONS)
    if not conv.ReadFile(mol, args_dict['inFileName']):
        raise ValueError("Error reading file " + args_dict['inFileName'] + " using OpenBabel")
    return mol

def read_ob_mol(obmol, args_dict):
    """
    :param obmol: OBmol molecule
    :param args_dict: dictionary of processed command line arguments
    :return: A list of Atoms
    """

    num_atoms = obmol.NumAtoms()
    atoms = []

    for i in range(num_atoms):
        obatom = obmol.GetAtom(i + 1)
		# get symbol by atomic number
        """ C++ code uses here an elements table:
        if (atom->GetAtomicNum() <= ELEMENTS.size() && atom->GetAtomicNum() > 0) {
			mol->_symbol[i] = strdup(ELEMENTS[atom->GetAtomicNum() - 1].c_str());
		}
		else {
			mol->_symbol[i] = strdup(atom->GetType());
		}"""

        # TODO: check if elements table is needed

        if args_dict["ignoreSym"]:
            symbol = "XX"
        else:
            symbol = obatom.GetType()
        position = (obatom.GetX(), obatom.GetY(), obatom.GetZ())

        atom = Atom(symbol, position, args_dict["useMass"])

        adjacent = []
        iter = openbabel.OBAtomAtomIter(obatom)
        for neighbour_atom in iter:
            adjacent.append(neighbour_atom.GetIdx())
        atom.adjacent = adjacent

        atoms.append(atom)
    #TODO: mol->initSimilarity(DEPTH_ITERATIONS);???
    return atoms


def read_csm_file(f, args_dict):
    """
    :param f: CSM file (the file object, not the file name)
    :param args_dict: dictionary of processed command line arguments
    :return: A list of Atoms
    """

    try:
        size = int(f.readline())
    except ValueError:
        raise ValueError("Input Error: Number of atoms not supplied")

    if size > 0:
        atoms = []
    else:
        return None

    for i in range(size):
        line = f.readline().split()
        try:
            if args_dict["ignoreSym"]:
                symbol = "XX"
            else:
                symbol = line[0]
            position = (float(line[1]), float(line[2]), float(line[3]))
            atom = Atom(symbol, position, args_dict["useMass"])
        except (ValueError, IndexError):
            raise ValueError("Input Error: Failed reading input for atom " + str(i+1))
        atoms.append(atom)

    for i in range(size):
        line = f.readline().split()
        try:
            atom_num = int(line.pop(0))
        except (ValueError, IndexError):
            raise ValueError("Input Error: Failed reading connectivity for atom " + str(i+1))
        if atom_num != i+1:
            raise ValueError("Input Error: Failed reading connectivity for atom " + str(i+1))

        neighbours = []
        for neighbour_str in line:
            try:
                neighbour = int(neighbour_str)
            except ValueError:
                raise ValueError("Input Error: Failed reading input for atom " + str(i+1))
            if neighbour > size:
                raise ValueError("Input Error: Failed reading input for atom " + str(i+1))
            neighbours.append(neighbour)

        atoms[i].adjacent = neighbours

    return atoms


def read_dir_file(f):
    """
    Reads a symmetry direction file
    :param f: File object
    :return: (x,y,z) of the symmetry axis
    """
    line = f.readline().split()
    result = (float(line[0]), float(line[1]), float(line[2]))
    return result


def read_perm_file(f):
    """
    Reads a permutation
    :param f: File object
    :return: permutation as a list of numbers

    Check that the permutation is legal, raise ValueError if not
    """
    line = f.readline().split()
    used = []
    for i in range(len(line)):
        used.append(False)

    result = []
    for num_str in line:
        try:
            num = int(num_str)
        except ValueError:
            raise ValueError("Invalid permutation")
        if num < 1 or num > len(line) or used[num-1]:
            raise ValueError("Invalid permutation")
        result.append(num-1)
        used[num-1] = True
    return result


if __name__ == '__main__':
    print("Testing the file reading functions")
    try:
        # file = open("../../test_cases/test1/AgCu10p1.xyz", "r")
        # file = open("../../test_cases/test3/c_in_1282148276_benzene.mol", "r")
        # file = open("../../test_cases/test6/AuCu10p6.xyz", "r")
        args_dict = {"useformat": False, "babelBond": False,
                     "inFileName": "../../test_cases/test3/c_in_1282148276_benzene.mol",
                     "ignoreSym": False, "useMass": True}
        obmol = open_non_csm_file(args_dict)
        result = read_ob_mol(obmol, args_dict)
        for a in result:
            print(a)
        # file.close()
    except IOError:
        print("no file")