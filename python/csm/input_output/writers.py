__author__ = 'YAEL'

import openbabel
from openbabel import OBConversion
import math
from csm.calculations.basic_calculations import check_perm_structure

def non_negative_zero(number):
    if math.fabs(number)<0.00001:
        return 0.0000
    else:
        return number

def print_results(result, in_args, calc_args, out_args):
    """
    Prints the CSM calculation results
    :param result: The result of the CSM calculation (a CSMState)
    :param in_args: Input arguments to CSM
    :param calc_args: Calculation arguments to CSM
    :param out_args: Output arguments to CSM
    """

    if calc_args['calc_type']=='just_perms':
        print("NUMBER OF PERMUTATIONS: %.4g" % result)
        return
    with open(out_args['out_file_name'], 'w', encoding='utf-8') as f:
        f.write("%s: %.4lf\n" % (calc_args['op_name'], abs(result.csm)))
        f.write("SCALING FACTOR: %7lf\n" % non_negative_zero(result.d_min))

        # print CSM, initial molecule, resulting structure and direction according to format specified
        try:
            percent_structure = check_perm_structure(result.molecule, result.perm)
            print("The permutation found maintains",
              str(round(percent_structure * 100, 2)) + "% of the original molecule's structure\n")
        except ValueError:
            print("The input molecule does not have bond information and therefore conservation of structure cannot be measured")
        if in_args['format'].lower() == "csm":
            print_output_csm(f, result, calc_args)
        else:
            print_output_ob(f, result, in_args, calc_args, out_args)



        # print norm

        if out_args['print_norm']:
            print("NORMALIZATION FACTOR: %7lf" % non_negative_zero(result.molecule.norm_factor))
            print("SCALING FACTOR OF SYMMETRIC STRUCTURE: %7lf" % non_negative_zero(result.d_min))
            print("DIRECTIONAL COSINES: %lf %lf %lf" % (non_negative_zero(result.dir[0]),
                                                        non_negative_zero(result.dir[1]),
                                                        non_negative_zero(result.dir[2])))
            print("NUMBER OF EQUIVALENCE GROUPS: %d" % len(result.molecule.equivalence_classes))

        # print local CSM

        if out_args['print_local']:
            sum = 0
            f.write("\nLocal CSM: \n")
            size = len(result.molecule.atoms)
            for i in range(size):
                sum += result.local_csm[i]
                f.write("%s %7lf\n" % (result.molecule.atoms[i].symbol, non_negative_zero(result.local_csm[i])))
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


def print_output_csm(f, result, calc_args):
    """
    Prints output in CSM format
    :param f: File to print to
    :param result: The result of the CSM calculation (a CSMState)
    :param calc_args: Calculation arguments to CSM
    """
    size = len(result.molecule.atoms)

    # print initial molecule

    f.write("\nINITIAL STRUCTURE COORDINATES\n%i\n\n" % size)
    for i in range(size):
        f.write("%3s%10lf %10lf %10lf\n" %
                (result.molecule.atoms[i].symbol,
                 non_negative_zero(result.molecule.atoms[i].pos[0]),
                 non_negative_zero(result.molecule.atoms[i].pos[1]),
                 non_negative_zero(result.molecule.atoms[i].pos[2])))

    for i in range(size):
        f.write("%d " % (i + 1))
        for j in result.molecule.atoms[i].adjacent:
            f.write("%d " % (j + 1))
        f.write("\n")

    # print resulting structure coordinates

    f.write("\nMODEL 02 RESULTING STRUCTURE COORDINATES\n%i\n" % size)

    for i in range(size):
        f.write("%3s%10lf %10lf %10lf\n" %
                (result.molecule.atoms[i].symbol,
                 non_negative_zero(result.symmetric_structure[i][0]),
                 non_negative_zero(result.symmetric_structure[i][1]),
                 non_negative_zero(result.symmetric_structure[i][2])))

    for i in range(size):
        f.write("%d " % (i + 1))
        for j in result.molecule.atoms[i].adjacent:
            f.write("%d " % (j + 1))
        f.write("\n")

    # print dir

    f.write("\n DIRECTIONAL COSINES:\n\n")
    f.write("%lf %lf %lf\n" % (non_negative_zero(result.dir[0]), non_negative_zero(result.dir[1]), non_negative_zero(result.dir[2])))

    print("%s: %.6lf" % (calc_args['op_name'], abs(result.csm)))
    print("CSM by formula: %.6lf" % (result.formula_csm))


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
    if str.lower(in_args['format'])=='pdb':
        f.write("\nMODEL 01")
    f.write("\nINITIAL STRUCTURE COORDINATES\n")

    num_atoms = result.molecule.obmol.NumAtoms()
    # update coordinates
    for i in range(num_atoms):
        try:
            atom = result.molecule.obmol.GetAtom(i + 1)
            atom.SetVector(non_negative_zero(result.molecule.atoms[i].pos[0]),
                       non_negative_zero(result.molecule.atoms[i].pos[1]),
                       non_negative_zero(result.molecule.atoms[i].pos[2]))
        except:
            pass

    write_ob_molecule(result.molecule.obmol, in_args['format'], f)

    # print resulting structure coordinates

    # update coordinates
    mol = result.molecule.obmol
    for i in range(num_atoms):
        try:
            a=mol.GetAtom(i+1)
            a.SetVector(non_negative_zero(result.symmetric_structure[i][0]),
                       non_negative_zero(result.symmetric_structure[i][1]),
                       non_negative_zero(result.symmetric_structure[i][2]))
        except:
            pass

    if str.lower(in_args['format'])=='pdb':
        f.write("\nMODEL 02")
    f.write("\nRESULTING STRUCTURE COORDINATES\n")
    write_ob_molecule(mol, in_args['format'], f)
    if str.lower(in_args['format']) == 'pdb':
        f.write("END\n")

    # print dir

    f.write("\n DIRECTIONAL COSINES:\n\n")
    f.write("%lf %lf %lf\n" % (non_negative_zero(result.dir[0]),
                               non_negative_zero(result.dir[1]),
                               non_negative_zero(result.dir[2])))

    if out_args['write_openu']:
        print("SV* %.6lf *SV" % abs(result.csm))
    else:
        print("%s: %.6lf" % (calc_args['op_name'], abs(result.csm)))
    print("CSM by formula: %.6lf" % (result.formula_csm))



def update_coordinates(obmol, out_atoms):
    """
    Updates the coordinates of the OpenBabel Molecule according to the Molecule data
    :param obmol: The OpenBable molecule
    :param outAtoms: The output atoms' coordinates
    """
    num_atoms = obmol.NumAtoms()
    for i in range(num_atoms):
        atom = obmol.GetAtom(i + 1)
        atom.SetVector(non_negative_zero(out_atoms[i]['pos'][0]), non_negative_zero(out_atoms[i]['pos'][1]), non_negative_zero(out_atoms[i]['pos'][2]))


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
    if str.lower(format)=='pdb':
        s=s.replace("END", "ENDMDL")
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

