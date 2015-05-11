from input_output.molecule import Atom

__author__ = 'YAEL'

def read_csm_file(f):
    """
    :param f: CSM file (the file object, not the file name)
    :return: A list of Atoms
    """

    try:
        size = int(f.readline())
    except ValueError:
        raise ValueError("Input Error: Number of atoms not supplied")

    if size > 0:
        atoms = []
        connectivity = []
    else:
        return None

    for i in range(size):
        line = f.readline().split()
        try:
            atom = (line[0], (float(line[1]), float(line[2]), float(line[3])))
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

        connectivity.append(neighbours)

    return atoms, connectivity


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
    file = open("../../test_cases/ALA.csm", "r")
    print(read_perm_file(file))
    file.close()