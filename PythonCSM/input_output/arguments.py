"""
Parse the CSM command line arguments.
"""
from argparse import ArgumentParser
import csv
from old_input_output.readers import read_dir_file, read_perm_file, read_csm_file, open_non_csm_file, read_ob_mol
from molecule.molecule import Molecule

__author__ = 'zmbq'


def _create_parser():
    parser = ArgumentParser()

    # The first three positional arguments
    parser.add_argument('type',
                        choices=('c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c8', 's2', 's4', 's6', 's8', 'cs', 'ci', 'ch'),
                        help='The type of operation')
    parser.add_argument('input', help='Input file')
    parser.add_argument('output', help='Output file')

    # Optional arguments (their names start with --)
    parser.add_argument('--ignoreHy', action='store_true', default=False, help='Ignore Hydrogen atoms in computations')
    parser.add_argument('--removeHy', action='store_true', default=False,
                        help='Remove Hydrogen atoms in computations, rebuild molecule without them and compute')
    parser.add_argument('--ignoreSym', action='store_true', default=False,
                        help='Ignore all atomic symbols, performing a purely geometric operation')

    parser.add_argument('--format', help='Use a specific input/output format')
    parser.add_argument('--writeOpenu', action='store_true', default=False,
                        help='Write output in open university format')
    parser.add_argument('--nolimit', action='store_true', default=False,
                        help='Allows running program while ignoring computational complexity')
    parser.add_argument('--useperm', type=str, help='Only compute for a single permutation')
    parser.add_argument('--usedir', type=str, help='Use a predefined axis as a starting point. '
                                                   'This options ignores the -ignoreSym/-ignoreHy/-removeHy flags')
    parser.add_argument('--findperm', action='store_true', default=False, help='Attempt to search for a permutation')
    parser.add_argument('--detectOutliers', action='store_true', default=False,
                        help="Use statistical methods to try and improve --findperm's results")
    parser.add_argument('--babelbond', action='store_true', default=False, help='Let OpenBabel compute bonding')
    parser.add_argument('--useMass', action='store_true', default=False,
                        help='Use the atomic masses to define center of mass')
    parser.add_argument('--timeOnly', action='store_true', default=False, help="Only print the time and exit")
    parser.add_argument('--babelTest', action='store_true', default=False, help="Test if the molecule is legal or not")
    parser.add_argument('--sn_max', type=int, help='The maximal sn to try, relevant only for chirality')
    parser.add_argument('--printNorm', action='store_true', default=False, help='Print the normalization factor as well')
    parser.add_argument('--printLocal', action='store_true', default=False,
                        help='Print the local CSM (csm for each atom) in the output file')
    parser.add_argument('--approx', action='store_true', default=False,
                        help='Equivalent to --detectOutliers --findperm together')
    parser.add_argument('--keepCenter', action='store_true', default=False,
                        help='Do not change coordinates s.t. (0,0,0) corresponds to Center of Mass')
    parser.add_argument('--log', type=str, help='Write a detailed log to logfile')
    parser.add_argument('--outputPerms', action='store',
                        help='Writes all enumerated permutations to file')
    parser.add_argument('--useChains', action='store_true', default=False,
                        help='Use chains specified in the PDB file in order to calculate permutations')

    return parser


def _check_arguments(processed):

    if processed['sn_max'] and processed['type'] != 'CH':
        raise ValueError("Option --sn_max only applies to chirality")

    if (processed['findPerm'] and 'perm' in processed) or (processed['findPerm'] and 'dirFile' in processed) \
            or ('dirFile' in processed and 'perm' in processed):
        raise ValueError("--findperm, --useperm and --usedir are mutually exclusive")

    if processed['removeHy'] and processed['ignoreHy']:
        raise ValueError("--removeHy and --ignoreHy are mutually exclusive")

    if "perm" in processed:

        if processed['type'] == 'CH':
            raise ValueError("Chirality can't be given a permutation, run the specific csm operation instead")

        # In C++ code ignoreSym, ignoreHy and removeHy are used only when usePerm is false
        if processed["ignoreSym"]:
            raise ValueError("--useperm ignores the --ignoreSym option, can't use them together")
        if processed["ignoreHy"]:
            raise ValueError("--useperm ignores the --ignoreHy option, can't use them together")
        if processed["removeHy"]:
            raise ValueError("--useperm ignores the --removeHy option, can't use them together")

        if len(processed["perm"]) != len(processed['molecule'].atoms):
            raise ValueError("Invalid permutation")

    if processed['useChains'] and not processed['molecule'].chains:
        raise ValueError("--useChains specified but no chains provided in the molecule file")


def _open_files(parse_res, result):

    # try to open and read the input file
    try:
        if result['format'].lower() == "csm":
            with open(parse_res.input, 'r') as infile:
                atoms = read_csm_file(infile, result)
                result['molecule'] = Molecule(atoms)
        else:
            result["obmol"] = open_non_csm_file(result)
            (atoms, chains) = read_ob_mol(result["obmol"], result)
            if not parse_res.useChains:  # Todo: Print chains
                chains = []
            result['molecule'] = Molecule(atoms, chains=chains, obmol=result["obmol"])
    except IOError:
        raise ValueError("Failed to open data file " + parse_res.input)
    except ValueError:
        raise ValueError("Failed to read molecule from data file")

    # try to open the output file for writing
    try:
        result['outFile'] = open(parse_res.output, 'w')
    except IOError:
        raise ValueError("Failed to open output file " + parse_res.output + " for writing")

    # try to open the permFile for reading (if exists)
    if parse_res.useperm:
        try:
            with open(parse_res.useperm, 'r') as permfile:
                perm = read_perm_file(permfile)
                result['perm'] = perm
        except IOError:
            raise ValueError("Failed to open perm file " + parse_res.useperm + " for reading")

    # try to open the dirFile for reading (if exists)
    if parse_res.usedir:
        try:
            with open(parse_res.usedir, 'r') as dirfile:
                dir = read_dir_file(dirfile)
                result['dir'] = dir
        except IOError:
            raise ValueError("Failed to open dir file " + parse_res.usedir + " for reading")
        except (ValueError, IndexError):
            raise ValueError("Can't read legal direction from file " + parse_res.usedir)

    # Try to open the permutation output file
    if parse_res.outputPerms:
        try:
            result['outPermFile'] = open(parse_res.outputPerms, 'w')
        except IOError:
            raise ValueError("Failed to open output permutation file " + parse_res.outputPerms + " for writing")
    else:
        result['outPermFile'] = None





def _process_arguments(parse_res):

    result = {}

    op_names = {
        "cs": ('CS', 2, "MIRROR SYMMETRY"),
        "ci": ('CI', 2, "INVERSION (S2)"),
        "ch": ('CH', 2, "CHIRALITY"),
        "c2": ('CN', 2, "C2 SYMMETRY"),
        'c3': ('CN', 3, "C3 SYMMETRY"),
        'c4': ('CN', 4, "C4 SYMMETRY"),
        'c5': ('CN', 5, "C5 SYMMETRY"),
        'c6': ('CN', 6, "C6 SYMMETRY"),
        'c7': ('CN', 7, "C7 SYMMETRY"),
        'c8': ('CN', 8, "C8 SYMMETRY"),
        's2': ('SN', 2, "S2 SYMMETRY"),
        's4': ('SN', 4, "S4 SYMMETRY"),
        's6': ('SN', 6, "S6 SYMMETRY"),
        's8': ('SN', 8, "S8 SYMMETRY")
    }

    result['type'] = op_names[parse_res.type][0]
    result['opOrder'] = op_names[parse_res.type][1]
    result['opName'] = op_names[parse_res.type][2]
    result['inFileName'] = parse_res.input
    result['outFileName'] = parse_res.output

    result['sn_max'] = parse_res.sn_max

    result['ignoreHy'] = parse_res.ignoreHy
    result['removeHy'] = parse_res.removeHy
    result['ignoreSym'] = parse_res.ignoreSym

    result['format'] = parse_res.format
    result['useformat'] = result['format'] is not None
    if not result['format']:
        # get input file extension
        result['format'] = parse_res.input.split(".")[-1]

    result['writeOpenu'] = parse_res.writeOpenu
    result['limitRun'] = not parse_res.nolimit
    result['babelBond'] = parse_res.babelbond
    result['useMass'] = parse_res.useMass
    result['timeOnly'] = parse_res.timeOnly
    result['printNorm'] = parse_res.printNorm
    result['findPerm'] = parse_res.findperm
    result['detectOutliers'] = parse_res.detectOutliers
    result['useChains'] = parse_res.useChains
    if parse_res.approx:
        result['findPerm'] = True
        result['detectOutliers'] = True
    result['babelTest'] = parse_res.babelTest
    result['printLocal'] = parse_res.printLocal
    result['keepCenter'] = parse_res.keepCenter
    if parse_res.writeOpenu:
        result['format'] = "PDB"
    result['logFileName'] = parse_res.log

    #_open_files(parse_res, result)
    _check_arguments(result)

    if not result['sn_max']:
        result['sn_max'] = 8

    return result

def get_arguments(args):
    """
    Parse the command line arguments and return a dictionary with all the CSM run options
    :param args: Command line arguments (sys.args, for example)
    :return: A dictionary with all the CSM run options.
    """
    parser = _create_parser()
    parsed_args = parser.parse_args(args)
    csm_args = _process_arguments(parsed_args)

    return csm_args

def _process_split_arguments(parse_res):

    mol_result = {}
    calc_result= {}
    out_result= {}

    op_names = {
        "cs": ('CS', 2, "MIRROR SYMMETRY"),
        "ci": ('CI', 2, "INVERSION (S2)"),
        "ch": ('CH', 2, "CHIRALITY"),
        "c2": ('CN', 2, "C2 SYMMETRY"),
        'c3': ('CN', 3, "C3 SYMMETRY"),
        'c4': ('CN', 4, "C4 SYMMETRY"),
        'c5': ('CN', 5, "C5 SYMMETRY"),
        'c6': ('CN', 6, "C6 SYMMETRY"),
        'c7': ('CN', 7, "C7 SYMMETRY"),
        'c8': ('CN', 8, "C8 SYMMETRY"),
        's2': ('SN', 2, "S2 SYMMETRY"),
        's4': ('SN', 4, "S4 SYMMETRY"),
        's6': ('SN', 6, "S6 SYMMETRY"),
        's8': ('SN', 8, "S8 SYMMETRY")
    }

    calc_result['type'] = op_names[parse_res.type][0]
    calc_result['op_order'] = op_names[parse_res.type][1]
    calc_result['op_name'] = op_names[parse_res.type][2]
    calc_result['sn_max'] = parse_res.sn_max
    calc_result['limit_run'] = not parse_res.nolimit
    calc_result['find_perm'] = parse_res.findperm
    calc_result['detect_outliers'] = parse_res.detectOutliers
    if parse_res.approx:
        calc_result['approx']=True
        calc_result['find_perm'] = True
        calc_result['detect_outliers'] = True

    mol_result['in_file_name'] = parse_res.input
    mol_result['ignore_hy'] = parse_res.ignoreHy
    mol_result['remove_hy'] = parse_res.removeHy
    mol_result['ignore_sym'] = parse_res.ignoreSym
    mol_result['format'] = parse_res.format
    mol_result['useformat'] = mol_result['format'] is not None
    if not mol_result['format']:
        # get input file extension
        mol_result['format'] = parse_res.input.split(".")[-1]
    mol_result['babel_bond'] = parse_res.babelbond
    mol_result['use_mass'] = parse_res.useMass
    mol_result['use_chains'] = parse_res.useChains
    mol_result['keep_center'] = parse_res.keepCenter
    if parse_res.writeOpenu:
        mol_result['format'] = "PDB"

    out_result['write_openu'] = parse_res.writeOpenu
    out_result['print_norm'] = parse_res.printNorm
    out_result['print_local'] = parse_res.printLocal
    out_result['log_file_name'] = parse_res.log
    out_result['out_file_name'] = parse_res.output

    _open_files(parse_res, mol_result)
    _check_arguments(mol_result)
    _check_arguments(calc_result)
    _check_arguments(out_result)

    if not calc_result['sn_max']:
        calc_result['sn_max'] = 8

    return mol_result, calc_result, out_result

def get_split_arguments(args):
    """
    :param args:
    :return:
    """
    parser=_create_parser()
    parsed_args = parser.parse_args(args)
    mol_args, calc_args, out_args=_process_split_arguments(parsed_args)
    return mol_args, calc_args, out_args