import json

from csm.input_output.formatters import format_CSM
import io

__author__ = 'YAEL'

import openbabel
from openbabel import OBConversion
import math
from csm.calculations.basic_calculations import check_perm_structure, check_perm_cycles
from csm.molecule.molecule import MoleculeFactory

def non_negative_zero(number):
    if math.fabs(number)<0.00001:
        return 0.0000
    else:
        return number


def json_results(result, dictionary_args, result_string):
    json_dict={"Result":
        {
        "molecule": result.molecule.to_json(),
        "op_order": result.op_order,
        "op_type": result.op_type,
        "csm": result.csm,

        "result_string":result_string,
        "perm": result.perm,
        "dir": list(result.dir),
        "d_min": result.d_min,
        "symmetric_structure": [list(i) for i in result.symmetric_structure],
        "local_csm": result.local_csm,
        "perm_count": result.perm_count,
        "formula_csm": result.formula_csm,
        "normalized_molecule_coords": [list(i) for i in result.normalized_molecule_coords],
        "normalized_symmetric_structure": [list(i) for i in result.normalized_symmetric_structure],
             }
           }
    with open(dictionary_args['out_file_name'], 'w', encoding='utf-8') as f:
        json.dump(json_dict, f)


def _print_results(f, result, dictionary_args):
    f.write("%s: %.4lf\n" % (dictionary_args['op_name'], abs(result.csm)))
    f.write("SCALING FACTOR: %7lf\n" % non_negative_zero(result.d_min))

    # print CSM, initial molecule, resulting structure and direction according to format specified
    try:
        percent_structure = check_perm_structure(result.molecule, result.perm)
        print("The permutation found maintains",
              str(round(percent_structure * 100, 2)) + "% of the original molecule's structure")

    except ValueError:
        print(
            "The input molecule does not have bond information and therefore conservation of structure cannot be measured")

    falsecount, num_invalid, cycle_counts = check_perm_cycles(result.perm, dictionary_args['op_order'],
                                                              dictionary_args['op_type'])
    if falsecount > 0 or dictionary_args['calc_type'] == 'approx':
        print("The permutation found contains %d invalid %s. %.2lf%% of the molecule's atoms are in legal cycles" % (
        falsecount, "cycle" if falsecount == 1 else "cycles",
        100 * (len(result.molecule) - num_invalid) / len(result.molecule)))
        for cycle_len in sorted(cycle_counts):
            valid = cycle_len == 1 or cycle_len == dictionary_args['op_order'] or (
            cycle_len == 2 and dictionary_args['op_type'] == 'SN')
            count = cycle_counts[cycle_len]
            print("There %s %d %s %s of length %d" % (
            "is" if count == 1 else "are", count, "invalid" if not valid else "", "cycle" if count == 1 else "cycles",
            cycle_len))

    if dictionary_args['format'].lower() == "csm":
        print_output_csm(f, result, dictionary_args)
    else:
        print_output_ob(f, result, dictionary_args, dictionary_args, dictionary_args)

    # print local CSM

    if dictionary_args['print_local']:
        sum = 0
        f.write("\nLocal CSM: \n")
        size = len(result.molecule.atoms)
        for i in range(size):
            sum += result.local_csm[i]
            f.write("%s %7lf\n" % (result.molecule.atoms[i].symbol, non_negative_zero(result.local_csm[i])))
        f.write("\nsum: %7lf\n" % sum)

    # print chirality
    if dictionary_args['op_type'] == 'CH':
        if result.op_type == 'CS':
            f.write("\n MINIMUM CHIRALITY WAS FOUND IN CS\n\n")
        else:
            f.write("\n MINIMUM CHIRALITY WAS FOUND IN S%d\n\n" % result.op_order)

    # print permutation
    f.write("\n PERMUTATION:\n\n")
    for i in result.perm:
        f.write("%d " % (i + 1))
    f.write("\n")

def print_results(result, dictionary_args):
    """
    Prints the CSM calculation results
    :param result: The result of the CSM calculation (a CSMState)
    :param dictionary_args: Input arguments to CSM
    :param dictionary_args: Calculation arguments to CSM
    :param dictionary_args: Output arguments to CSM
    """


    if dictionary_args['calc_type']== 'just_perms':
        print("NUMBER OF PERMUTATIONS: %5.4g" % result)
        return

    if dictionary_args['json_output']:
        result_io=io.StringIO()
        _print_results(result_io, result, dictionary_args)
        result_string=result_io.getvalue()
        result_io.close()
        json_results(result, dictionary_args, result_string)
        return



    with open(dictionary_args['out_file_name'], 'w', encoding='utf-8') as f:
        _print_results(f, result, dictionary_args)



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
        f.write("%3s%10.5lf %10.5lf %10.5lf\n" %
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
        f.write("%3s%10.5lf %10.5lf %10.5lf\n" %
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

    print("%s: %s" % (calc_args['op_name'], format_CSM(result.csm)))
    print("CSM by formula: %s" % format_CSM(result.formula_csm))


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

    obmol=MoleculeFactory._obm_from_file(result.molecule._filename, result.molecule._format, result.molecule._babel_bond)[0]
    for to_remove in result.molecule._deleted_atom_indices:
        obmol.DeleteAtom(obmol.GetAtom(to_remove + 1))

    num_atoms = obmol.NumAtoms()
    # update coordinates
    for i in range(num_atoms):
        try:
            atom = obmol.GetAtom(i + 1)
            atom.SetVector(non_negative_zero(result.molecule.atoms[i].pos[0]),
                       non_negative_zero(result.molecule.atoms[i].pos[1]),
                       non_negative_zero(result.molecule.atoms[i].pos[2]))
        except Exception as e:
            pass

    write_ob_molecule(obmol, in_args['format'], f)

    # print resulting structure coordinates

    # update output coordinates to match symmetric structure
    for i in range(num_atoms):
        try:
            a=obmol.GetAtom(i+1)
            a.SetVector(non_negative_zero(result.symmetric_structure[i][0]),
                       non_negative_zero(result.symmetric_structure[i][1]),
                       non_negative_zero(result.symmetric_structure[i][2]))
        except Exception as e:
            pass

    if str.lower(in_args['format'])=='pdb':
        f.write("\nMODEL 02")
    f.write("\nRESULTING STRUCTURE COORDINATES\n")
    write_ob_molecule(obmol, in_args['format'], f)
    if str.lower(in_args['format']) == 'pdb':
        f.write("END\n")

    # print dir

    f.write("\n DIRECTIONAL COSINES:\n\n")
    f.write("%lf %lf %lf\n" % (non_negative_zero(result.dir[0]),
                               non_negative_zero(result.dir[1]),
                               non_negative_zero(result.dir[2])))

    # if out_args['write_openu']:
    #     print("SV* %.6lf *SV" % abs(result.csm))
    # else:
    print("%s: %.6lf" % (calc_args['op_name'], abs(result.csm)))
    print("CSM by formula: %.6lf" % (result.formula_csm))


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
