"""
Parse the CSM command line arguments.
"""
from argparse import ArgumentParser

import sys

from collections import namedtuple

__author__ = 'zmbq'


class OurParser(ArgumentParser):
    def error(self, message):
        print("Error: %s" % message, file=sys.stderr)
        print("Enter csm --help for help", file=sys.stderr)
        sys.exit(2)


def _create_parser():
    parser = OurParser(usage="\ncsm type input_molecule output_file [additional arguments]")

    # The first three positional arguments
    c_symmetries = ['c%d' % n for n in range(2, 21)]
    s_symmetries = ['s%d' % n for n in range(2, 21, 2)]
    parser.add_argument('type',
                        choices=c_symmetries + s_symmetries + ['cs', 'ci', 'ch'],
                        help='The type of operation')
    parser.add_argument('input', help='Input molecule file')
    parser.add_argument('output', default='output.txt', help='Output file')

    # Optional arguments (their names start with --)
    parser.add_argument('--trivial', action='store_true', default=False,
                        help='CSM of identity perm, or, if chains, CSM of chain permutation with no atom permutation')
    parser.add_argument('--just-perms', action='store_true', default=False,
                        help='no calculation of CSM. without --outputPerms, only counts the perm. ')
    parser.add_argument('--approx', action='store_true', default=False,
                        help='Equivalent to --detectOutliers --findperm together')

    parser.add_argument('--use-perm', type=str, help='Only compute for a single permutation')
    parser.add_argument('--use-dir', type=str, help='Use a predefined axis as a starting point. '
                                                   'This options ignores the --ignore-sym/--ignore-hy/--remove-hy flags')
    parser.add_argument('--find-perm', action='store_true', default=False, help='Attempt to search for a permutation')
    parser.add_argument('--detect-outliers', action='store_true', default=False,
                        help="Use statistical methods to try and improve --findperm's results")

    parser.add_argument('--keep-structure', action='store_true', default=False,
                        help='Maintain molecule structure from being distorted')
    parser.add_argument('--ignore-hy', action='store_true', default=False, help='Ignore Hydrogen atoms in computations')
    parser.add_argument('--remove-hy', action='store_true', default=False,
                        help='Remove Hydrogen atoms in computations, rebuild molecule without them and compute')
    parser.add_argument('--ignore-sym', action='store_true', default=False,
                        help='Ignore all atomic symbols, performing a purely geometric operation')

    parser.add_argument('--format', help='Use a specific input/output format')
    parser.add_argument('--write-openu', action='store_true', default=False,
                        help='Write output in open university format')
    parser.add_argument('--no-limit', action='store_true', default=False,
                        help='Allows running program while ignoring computational complexity')

    parser.add_argument('--babel-bond', action='store_true', default=False, help='Let OpenBabel compute bonding')
    parser.add_argument('--no-babel',  action='store_true', default=False, help='force suppress automatically using babelbond to create bonds')
    parser.add_argument('--use-mass', action='store_true', default=False,
                        help='Use the atomic masses to define center of mass')
    parser.add_argument('--time-only', action='store_true', default=False, help="Only print the time and exit")
    parser.add_argument('--babel-test', action='store_true', default=False, help="Test if the molecule is legal or not")
    parser.add_argument('--sn-max', type=int, default=8, help='The maximal sn to try, relevant only for chirality')
    parser.add_argument('--print-norm', action='store_true', default=False,
                        help='Print the normalization factor as well')
    parser.add_argument('--print-local', action='store_true', default=False,
                        help='Print the local CSM (csm for each atom) in the output file')

    parser.add_argument('--log', type=str, help='Write a detailed log to logfile')
    parser.add_argument('--output-perms', action='store', default=None,
                        help='Writes all enumerated permutations to file')
    parser.add_argument('--use-chains', action='store_true', default=False,
                        help='Use chains specified in the PDB file in order to calculate permutations')
    parser.add_argument('--hungarian', action='store_true', default=False,
                        help='Use hungarian algorithm in approx')
    parser.add_argument('--print-approx', action='store_true', default=False,
                        help='add some printouts to approx')

    return parser


def _check_arguments(in_args, calc_args, out_args):
    if (calc_args['find_perm'] and 'perm' in calc_args) or (calc_args['find_perm'] and 'dir_file_name' in calc_args) \
            or ('dir_file_name' in calc_args and 'perm' in calc_args):
        raise ValueError("--findperm, --useperm and --usedir are mutually exclusive")

    if in_args['remove_hy'] and in_args['ignore_hy']:
        raise ValueError("--removeHy and --ignoreHy are mutually exclusive")

    if "perm" in calc_args:
        if calc_args['type'] == 'CH':
            raise ValueError("Chirality can't be given a permutation, run the specific csm operation instead")

        # In C++ code ignoreSym, ignoreHy and removeHy are used only when usePerm is false
        if in_args["ignore_sym"]:
            raise ValueError("--useperm ignores the --ignoreSym option, can't use them together")
        if in_args["ignore_hy"]:
            raise ValueError("--useperm ignores the --ignoreHy option, can't use them together")
        if in_args["remove_hy"]:
            raise ValueError("--useperm ignores the --removeHy option, can't use them together")

    if calc_args['detect_outliers'] and not calc_args['find_perm']:
        raise ValueError("--detectOutliers must be used with --findperm")

    #if in_args['use_chains'] and not in_args['molecule'].chains:
     #   raise ValueError("--useChains specified but no chains provided in the molecule file")


OperationCode = namedtuple('OperationCode', ('type', 'order', 'name'))
_opcode_data = {
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
    'c10': ('CN', 10, "C10 SYMMETRY"),
    's2': ('SN', 2, "S2 SYMMETRY"),
    's4': ('SN', 4, "S4 SYMMETRY"),
    's6': ('SN', 6, "S6 SYMMETRY"),
    's8': ('SN', 8, "S8 SYMMETRY"),
    's10': ('SN', 8, "S10 SYMMETRY")
}


def get_operation_data(opcode):
    """
    Returns data about an operation based on the opcode
    Args:
        opcode: c2, s4, etc...

    Returns:
        And OperationCode object, with type, order and name
    """
    def isint(s):
        try:
            int(s)
            return True
        except ValueError:
            return False

    opcode = opcode.lower()
    if opcode[0]=='c' and isint(opcode[1:]):
        return OperationCode(type='CN', order=int(opcode[1:]), name=opcode.upper() + ' SYMMETRY')
    if opcode[0]=='s' and isint(opcode[1:]):
        return OperationCode(type='SN', order=int(opcode[1:]), name=opcode.upper() + ' SYMMETRY')
    try:
        data = _opcode_data[opcode.lower()]
    except KeyError:
        raise
    return OperationCode(type=data[0], order=data[1], name=data[2])


def _process_split_arguments(parse_res):
    """
    Divides the parsed arguments (from argparse) into three dictionaries
    Args:
        parse_res: Result of argparse

    Returns:
        (in_args, calc_args, out_args) - input arguments, calculation arguments and output arguments

    """

    mol_args = {}
    calc_args = {}
    out_args = {}

    op = get_operation_data(parse_res.type)

    calc_args['op_type'] = op.type
    calc_args['op_order'] = op.order
    calc_args['op_name'] = op.name
    calc_args['sn_max'] = parse_res.sn_max
    calc_args['limit_run'] = not parse_res.no_limit
    calc_args['print_approx']= parse_res.print_approx

    calc_args['keep_structure'] =    mol_args['keep_structure']= parse_res.keep_structure


    calc_args['calc_type'] = 'exact'  # this is the default, which will be changed if relevant

    #calc_args['just_perms'] = parse_res.justperms
    if parse_res.just_perms:
        calc_args['calc_type'] = 'just_perms'

    if parse_res.trivial:
        calc_args['calc_type'] = 'trivial'

    calc_args['detect_outliers'] = parse_res.detect_outliers
    if parse_res.find_perm:
        calc_args['calc_type'] = 'approx'
    # calc_args['find_perm'] = parse_res.findperm
    if parse_res.approx:
        calc_args['calc_type'] = 'approx'
        #calc_args['find_perm'] = True
        calc_args['detect_outliers'] = True
    if parse_res.hungarian:
        calc_args['hungarian'] = True

    if parse_res.output_perms:
        calc_args['print_perms'] = True

    mol_args['in_file_name'] = parse_res.input
    mol_args['ignore_hy'] = parse_res.ignore_hy
    mol_args['remove_hy'] = parse_res.remove_hy
    mol_args['ignore_symm'] = parse_res.ignore_sym
    mol_args['format'] = parse_res.format
    mol_args['useformat'] = mol_args['format'] is not None
    if not mol_args['format']:
        # get input file extension
        mol_args['format'] = parse_res.input.split(".")[-1]
    mol_args['babel_bond'] = parse_res.babel_bond
    mol_args['no_babel'] = parse_res.no_babel
    mol_args['use_mass'] = parse_res.use_mass
    calc_args['use_chains'] = mol_args['use_chains'] = parse_res.use_chains
    if parse_res.write_openu:
        mol_args['format'] = "PDB"
    if parse_res.use_perm:
        mol_args['perm_file_name'] = parse_res.use_perm
    if parse_res.use_dir:
        mol_args['dir_file_name'] = parse_res.use_dir

    out_args['write_openu'] = parse_res.write_openu
    out_args['print_norm'] = parse_res.print_norm
    out_args['print_local'] = parse_res.print_local
    out_args['log_file_name'] = parse_res.log
    out_args['out_file_name'] = parse_res.output
    out_args['perms_csv_name'] = parse_res.output_perms

    #_check_arguments(mol_args, calc_args, out_args)

    return mol_args, calc_args, out_args


def get_split_arguments(args):
    """
    :param args:
    :return:
    """
    parser = _create_parser()
    parsed_args = parser.parse_args(args)
    mol_args, calc_args, out_args = _process_split_arguments(parsed_args)
    return mol_args, calc_args, out_args
