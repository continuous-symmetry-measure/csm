__author__ = 'YAEL'

from calculations.molecule import Atom
from openbabel import OBConversion


def print_all_output(output_obj, args_dict):
    """
    Prints all the outputs
    :param output_obj: Object with all old_calculations outputs
    :param args_dict: Dictionary with all processed arguments
    """
    f = args_dict['outFile']
    f.write("%s: %.4lf\n" % (args_dict['opName'], abs(output_obj.csm)))
    f.write("SCALING FACTOR: %7lf\n" % output_obj.dMin)

    # print CSM, initial molecule, resulting structure and direction according to format specified

    if args_dict['format'].lower() == "csm":
        print_output(output_obj, args_dict)
    else:
        print_output_format(output_obj, args_dict)

    # print norm

    if args_dict['printNorm']:
        print("NORMALIZATION FACTOR: %7lf" % output_obj.molecule.norm_factor)
        print("SCALING FACTOR OF SYMMETRIC STRUCTURE: %7lf" % output_obj.dMin)
        print("DIRECTIONAL COSINES: %lf %lf %lf" % (output_obj.dir[0], output_obj.dir[1], output_obj.dir[2]))
        print("NUMBER OF EQUIVALENCE GROUPS: %d" % len(output_obj.molecule.equivalence_classes))

    # print local CSM

    if args_dict['printLocal']:
        sum = 0
        f.write("\nLocal CSM: \n")
        size = len(output_obj.molecule.atoms)
        for i in range(size):
            sum += output_obj.localCSM[i]
            f.write("%s %7lf\n" % (output_obj.molecule.atoms[i].symbol, output_obj.localCSM[i]))
        f.write("\nsum: %7lf\n" % sum)

    # print chirality

    if args_dict['type'] == 'CH':
        if output_obj.chMinType == 'CS':
            f.write("\n MINIMUM CHIRALITY WAS FOUND IN CS\n\n")
        else:
            f.write("\n MINIMUM CHIRALITY WAS FOUND IN S%d\n\n" % output_obj.chMinOrder)

    # print permutation

    f.write("\n PERMUTATION:\n\n")
    for i in output_obj.perm:
        f.write("%d " % (i + 1))
    f.write("\n")

    f.close()


def print_output(output_obj, args_dict):
    """
    Prints output in CSM format
    :param output_obj: Object with all old_calculations outputs
    :param args_dict: Dictionary with all processed arguments
    """

    f = args_dict['outFile']
    print("%s: %.6lf" % (args_dict['opName'], abs(output_obj.csm)))
    size = len(output_obj.molecule.atoms)

    # print initial molecule

    f.write("\n INITIAL STRUCTURE COORDINATES\n%i\n\n" % size)
    for i in range(size):
        f.write("%3s%10lf %10lf %10lf\n" %
                (output_obj.molecule.atoms[i].symbol,
                 output_obj.molecule.atoms[i].pos[0],
                 output_obj.molecule.atoms[i].pos[1],
                 output_obj.molecule.atoms[i].pos[2]))

    for i in range(size):
        f.write("%d " % (i + 1))
        for j in output_obj.molecule.atoms[i].adjacent:
            f.write("%d " % (j + 1))
        f.write("\n")

    # print resulting structure coordinates

    f.write("\n RESULTING STRUCTURE COORDINATES\n%i\n" % size)

    for i in range(size):
        f.write("%3s%10lf %10lf %10lf\n" %
                (output_obj.molecule.atoms[i].symbol,
                 output_obj.outAtoms[i][0],
                 output_obj.outAtoms[i][1],
                 output_obj.outAtoms[i][2]))

    for i in range(size):
        f.write("%d " % (i + 1))
        for j in output_obj.molecule.atoms[i].adjacent:
            f.write("%d " % (j + 1))
        f.write("\n")

    # print dir

    f.write("\n DIRECTIONAL COSINES:\n\n")
    f.write("%lf %lf %lf\n" % (output_obj.dir[0], output_obj.dir[1], output_obj.dir[2]))


def print_output_format(output_obj, args_dict):
    """
    Prints output using Open Babel
    :param output_obj: Object with all old_calculations outputs
    :param args_dict: Dictionary with all processed arguments
    """
    f = args_dict['outFile']

    # comment from C++:
    # TODO - should we print the centered molecule, or the original one (and, accordingly, the symmetric struct)

    # print initial molecule

    f.write("\n INITIAL STRUCTURE COORDINATES\n")

    num_atoms = args_dict['obmol'].NumAtoms()
    # update coordinates
    for i in range(num_atoms):
        atom = args_dict['obmol'].GetAtom(i + 1)
        atom.SetVector(output_obj.molecule.atoms[i].pos[0],
                       output_obj.molecule.atoms[i].pos[1],
                       output_obj.molecule.atoms[i].pos[2])

    write_ob_molecule(args_dict['obmol'], args_dict['format'], f)

    # print resulting structure coordinates

    # update coordinates
    for i in range(num_atoms):
        atom = args_dict['obmol'].GetAtom(i + 1)
        atom.SetVector(output_obj.outAtoms[i][0],
                       output_obj.outAtoms[i][1],
                       output_obj.outAtoms[i][2])

    f.write("\n RESULTING STRUCTURE COORDINATES\n")
    write_ob_molecule(args_dict['obmol'], args_dict['format'], f)

    # print dir

    f.write("\n DIRECTIONAL COSINES:\n\n")
    f.write("%lf %lf %lf\n" % (output_obj.dir[0], output_obj.dir[1], output_obj.dir[2]))

    if args_dict['writeOpenu']:
        print("SV* %.4lf *SV\n" % abs(output_obj.csm))
    else:
        print("%s: %.4lf\n" % (args_dict['opName'], abs(output_obj.csm)))



def update_coordinates(obmol, out_atoms):
    """
    Updates the coordinates of the OpenBabel Molecule according to the Molecule data
    :param obmol: The OpenBable molecule
    :param outAtoms: The output atoms' coordinates
    """
    num_atoms = obmol.NumAtoms()
    for i in range(num_atoms):
        atom = obmol.GetAtom(i + 1)
        atom.SetVector(out_atoms[i]['pos'][0], out_atoms[i]['pos'][1], out_atoms[i]['pos'][2])


def write_ob_molecule(mol, format, f):
    """
    Write an Open Babel molecule to file
    :param mol: The molecule
    :param format: The output format
    :param f: The file to write output to
    :param f_name: The file's name (for extension-finding purpose)
    """
    conv = OBConversion()
    if not conv.SetOutFormat(format):
        raise ValueError("Error setting output format to " + format)

    # write to file

    try:
        s = conv.WriteString(mol)
    except (TypeError, ValueError, IOError):
        raise ValueError("Error writing data file using OpenBabel")
    f.write(s)


def print_molecule(molecule, f):
    f.write("The molecule:\n")
    for i in range(len(molecule)):
        f.write("%d %3s\n" % (i, molecule[i].symbol))
        f.write("\tconnections: ")
        for j in molecule[i].adjacent:
            f.write("%d, " %j )
        f.write("\n")


def print_equivalence_classes(groups, f):
    groups_num = len(groups)
    for i in range(groups_num):
        f.write("Group %d: " % i)
        for j in groups[i]:
            f.write("%d, " % j)
        f.write('\n')



if __name__ == '__main__':
    f = open("../testFiles/output.txt", "w")

    a1 = Atom("H", (2.0, 3.0, 5.225))
    a2 = Atom("H", (1.4, 3.0, 5.225))
    a3 = Atom("O", (2.0, 7.40, 5.225))

    a1.adjacent = [[2], [0, 2], []]
    a2.adjacent = [[2], [0, 2], [1]]
    a3.adjacent = [[], [2], []]

    atoms = [a1, a2, a3]
    outAtoms = [(2.0, 3.0, 5), (2.0, 3.0, 5), (2.0, 3.0, 5)]
    dir = [1.0, 2.0, 3.0]
    out = {'csm': 4.56, 'dMin': 3, 'atoms': atoms, 'outAtoms': outAtoms, 'dir': dir, 'groupNum': 2, 'norm': 2.3,
           'localCSM': [23.3, 2.4, -0.21]}
    print_output(out,
                 {'outFile': f, 'opName': 'C4 SYMMETRY', 'printNorm': True, 'printLocal': True})
