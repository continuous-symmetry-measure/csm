import json
import sys
from pathlib import Path

import os

from csm.calculations.basic_calculations import check_perm_structure_preservation, check_perm_equivalence, \
    check_perm_cycles
from csm.input_output.formatters import csm_log as print
from csm.molecule.molecule import MoleculeReader, Molecule, select_mols


def read_molecules(**kwargs):
    input_name = kwargs["in_file_name"]
    if os.path.isdir(input_name):
        kwargs.pop('in_file_name')
        mols = []
        files = []
        for directory, subdirectories, files1 in os.walk(input_name):
            files += files1
            break

        files.sort()

        supported_formats = [".mol", ".pdb", ".sdf", ".xyz", ".csm", ".cif", ".sd", ".mmcif"]

        have_warned = False
        i = 0
        mol_format = Path(files[i]).suffix
        try:
            while mol_format not in supported_formats:
                i += 1
                mol_format = Path(files[i]).suffix
        except:  # this is simply a check for a warning to user and not important enough to crash over
            pass

        for file_name in files:
            if file_name == "sym.txt":
                continue
            mol_file = os.path.join(input_name, file_name)
            mol_format_2 = Path(mol_file).suffix
            if mol_format_2 != mol_format and mol_format_2 in supported_formats and not have_warned:
                print("WARNING: you are reading multiple formats of files. Result printing may have errors")
                have_warned = True
            try:
                mol = MoleculeReader.from_file(mol_file, **kwargs)
                mols.append(mol)
            except Exception as e:
                if mol_format_2 in supported_formats:
                    print("failed to create a molecule from", file_name, str(e))
                # else:
                #    print(file_name, "was not read")

        mols, seld_mols = select_mols(mols, kwargs)

    elif not os.path.isfile(input_name):
        test = os.getcwd()
        raise ValueError("invalid file/folder name for molecule: ", input_name)
    else:  # file
        mols = MoleculeReader.multiple_from_file(**kwargs)
    sys.stderr.flush()

    if len(mols)<1:
        raise ValueError("Failed to read any molecules, calculation cannot proceed")
    return mols


def read_mols_from_std_in():
    raw_json = read_from_sys_std_in()
    less_raw_json = json.loads(raw_json)
    molecules = [Molecule.from_dict(json_dict) for json_dict in less_raw_json]
    return molecules


def read(**dictionary_args):
    mols = read_molecules(**dictionary_args)
    sys.stdout.write(json.dumps([mol.to_dict() for mol in mols], indent=4))
    return mols


def check_perm_validity(mol, perm, **kwargs):
    '''
    Checks whether a user-input permutation is valid-- has legal cycle lengths, only switches between equivalence classes,
    and maintains molecule bonds. Warns for each violation.
    :param mol: the molecule being permuted
    :param perm: the permutation
    '''
    falsecount, num_invalid, cycle_counts, bad_indices = check_perm_cycles(perm, kwargs['operation'])
    if falsecount > 0:
        print("Warning: Permutation does not maintain cycle structure")
    if not check_perm_equivalence(mol, perm):
        print("Warning: Permutation contains switches between non-equivalent atoms")
    try:
        if check_perm_structure_preservation(mol, perm) < 1:
            print("Warning: Permutation does not preserve molecule structure")
    except ValueError:  # molecule has no structure
        pass


def read_perm_file(molecule, filename):
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
            if num < 1 or num > len(molecule.fixed_indexes) or molecule.fixed_indexes[num - 1] is None:
                raise ValueError("Invalid permutation %s - out of range" % line)
            num -= 1
            if num in used:
                raise ValueError("Invalid permutation %s - number appears twice" % line)
            result.append(num)
            used.add(num)

    return result


def read_chain_perm_file(molecule, filename, **kwargs):
    perms=[]
    if len(molecule.chains)<4:
        print("Warning: the molecule has less than 4 chains, using --input-chain-perm is not recommended")
    with open(filename, 'r') as f:
        for line in f:
            perm=line.split()
            try:
                index_perm=[molecule.chains.name_to_index(name) for name in perm]
            except KeyError:
                raise ValueError("Permutation contains invalid chain name: " + line)
            if index_perm:
                falsecount, num_invalid, cycle_counts, bad_indices = check_perm_cycles(index_perm, kwargs['operation'])
                if falsecount > 0:
                    raise ValueError("Chain Permutation does not maintain cycle structure: "+line)

                for f_i, t_i in enumerate(index_perm):
                    if len(molecule.chains[f_i])!=len(molecule.chains[t_i]):
                        raise ValueError("ERROR: Cannot permute two chains with non matching lengths")

                if len(index_perm)<len(molecule.chains):
                    raise ValueError("Chain permutation is too short")

                perms.append(index_perm)
    return perms

def read_perm(molecule, perm_file_name=None, chain_perm_file_name=None, **kwargs):
    '''
    :param molecule: the molecule whose permutation is being read (necessary for validity checks)
    :param perm_file_name: the name of the file where the permutation is stored (file permutation is written in indexes starting from 1)
    :return: permutation- a list of permuted indexes (starting with 0) (can be None)
    '''
    return_perm = None
    if perm_file_name:
        perm = read_perm_file(molecule, perm_file_name)
        if len(perm) != len(molecule):
            raise ValueError("Invalid permutation - permutation is of size %d but molecule has %d atoms" %
                             (len(perm), len(molecule)))
        return_perm = [molecule.fixed_indexes[i] for i in perm]
        check_perm_validity(molecule, return_perm, **kwargs)
    if chain_perm_file_name:
        perms=read_chain_perm_file(molecule, chain_perm_file_name, **kwargs)
        if len(perms)<1:
            perms=None
        return_perm=perms
    return return_perm


def read_dir_file(filename):
    """
    THIS FUNCTION HAS BEEN RETIRED AND IS NOT IN USE AT THE MOMENT
    Reads a symmetry direction file
    :param filename: Name of dir file
    :return: (x,y,z) of the symmetry axis
    """
    dirs = []
    try:
        with open(filename, 'r') as f:
            for fileline in f:
                line = fileline.split()
                result = (float(line[0]), float(line[1]), float(line[2]))
                dirs.append(result)
        if len(dirs) == 0:
            raise ValueError("provided dir file is empty")
        return dirs
    except Exception as e:
        raise ValueError("Invalid input for use-dir")


def read_from_sys_std_in():
    '''
    a wrapper around calls to sys.stdin.read: we first check that there's something being piped, and if not, raise an error
    (otherwise the program hands forever)
    :return:
    '''
    if not sys.stdin.isatty():
        raw_json = sys.stdin.read()
        return raw_json
    else:
        raise ValueError("No input in sys.stdin")
