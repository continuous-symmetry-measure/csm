__author__ = 'YAEL'

from openbabel import OBConversion


def print_results(result, in_args, calc_args, out_args):
    """
    Prints the CSM calculation results
    :param result: The result of the CSM calculation (a CSMState)
    :param in_args: Input arguments to CSM
    :param calc_args: Calculation arguments to CSM
    :param out_args: Output arguments to CSM
    """
    with open(out_args['out_file_name'], 'w', encoding='utf-8') as f:
        f.write("%s: %.4lf\n" % (calc_args['op_name'], abs(result.csm)))
        f.write("SCALING FACTOR: %7lf\n" % result.d_min)

        # print CSM, initial molecule, resulting structure and direction according to format specified

        if in_args['format'].lower() == "csm":
            print_output(f, result, calc_args)
        else:
            print_output_ob(f, result, in_args, calc_args, out_args)

        # print norm

        if out_args['print_norm']:
            print("NORMALIZATION FACTOR: %7lf" % result.molecule.norm_factor)
            print("SCALING FACTOR OF SYMMETRIC STRUCTURE: %7lf" % result.d_min)
            print("DIRECTIONAL COSINES: %lf %lf %lf" % (result.dir[0], result.dir[1], result.dir[2]))
            print("NUMBER OF EQUIVALENCE GROUPS: %d" % len(result.molecule.equivalence_classes))

        # print local CSM

        if out_args['print_local']:
            sum = 0
            f.write("\nLocal CSM: \n")
            size = len(result.molecule.atoms)
            for i in range(size):
                sum += result.localCSM[i]
                f.write("%s %7lf\n" % (result.molecule.atoms[i].symbol, result.localCSM[i]))
            f.write("\nsum: %7lf\n" % sum)

        # print chirality
        if calc_args['op_type'] == 'CH':
            if result.op_type == 'CS':
                f.write("\n MINIMUM CHIRALITY WAS FOUND IN CS\n\n")
            else:
                f.write("\n MINIMUM CHIRALITY WAS FOUND IN S%d\n\n" % result.op_order)

        # print permutation
        f.write("\n PERMUTATION:\n\n")
        for i in result.perm:
            f.write("%d " % (i + 1))
        f.write("\n")


def print_output(f, result, calc_args):
    """
    Prints output in CSM format
    :param f: File to print to
    :param result: The result of the CSM calculation (a CSMState)
    :param calc_args: Calculation arguments to CSM
    """

    print("%s: %.6lf" % (calc_args['op_name'], abs(result.csm)))
    size = len(result.molecule.atoms)

    # print initial molecule

    f.write("\n INITIAL STRUCTURE COORDINATES\n%i\n\n" % size)
    for i in range(size):
        f.write("%3s%10lf %10lf %10lf\n" %
                (result.molecule.atoms[i].symbol,
                 result.molecule.atoms[i].pos[0],
                 result.molecule.atoms[i].pos[1],
                 result.molecule.atoms[i].pos[2]))

    for i in range(size):
        f.write("%d " % (i + 1))
        for j in result.molecule.atoms[i].adjacent:
            f.write("%d " % (j + 1))
        f.write("\n")

    # print resulting structure coordinates

    f.write("\n RESULTING STRUCTURE COORDINATES\n%i\n" % size)

    for i in range(size):
        f.write("%3s%10lf %10lf %10lf\n" %
                (result.molecule.atoms[i].symbol,
                 result.symmetric_structure[i][0],
                 result.symmetric_structure[i][1],
                 result.symmetric_structure[i][2]))

    for i in range(size):
        f.write("%d " % (i + 1))
        for j in result.molecule.atoms[i].adjacent:
            f.write("%d " % (j + 1))
        f.write("\n")

    # print dir

    f.write("\n DIRECTIONAL COSINES:\n\n")
    f.write("%lf %lf %lf\n" % (result.dir[0], result.dir[1], result.dir[2]))


def print_output_ob(f, result, in_args, calc_args, out_args):
    """
    Prints output using Open Babel
    :param f: File to write to
    :param result: The result of the CSM calculation (a CSMState)
    :param in_args: Input arguments to CSM
    :param calc_args: Calculation arguments to CSM
    :param out_args: Output arguments to CSM
    """
    # print initial molecule
    f.write("\n INITIAL STRUCTURE COORDINATES\n")

    num_atoms = result.molecule.obmol.NumAtoms()
    # update coordinates
    for i in range(num_atoms):
        try:
            atom = result.molecule.obmol.GetAtom(i + 1)
            atom.SetVector(result.molecule.atoms[i].pos[0],
                       result.molecule.atoms[i].pos[1],
                       result.molecule.atoms[i].pos[2])
        except:
            pass

    write_ob_molecule(result.molecule.obmol, in_args['format'], f)

    # print resulting structure coordinates

    # update coordinates
    for i in range(num_atoms):
        try:
            atom.SetVector(result.symmetric_structure[i][0],
                       result.symmetric_structure[i][1],
                       result.symmetric_structure[i][2])
        except:
            pass

    f.write("\n RESULTING STRUCTURE COORDINATES\n")
    write_ob_molecule(result.molecule.obmol, in_args['format'], f)

    # print dir

    f.write("\n DIRECTIONAL COSINES:\n\n")
    f.write("%lf %lf %lf\n" % (result.dir[0], result.dir[1], result.dir[2]))

    if out_args['write_openu']:
        print("SV* %.4lf *SV\n" % abs(result.csm))
    else:
        print("%s: %.4lf\n" % (calc_args['op_name'], abs(result.csm)))



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

