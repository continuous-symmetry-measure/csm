"""
Parse the CSM command line arguments.
"""
from argparse import ArgumentParser
import logging

from csm.calculations import permuters
from csm.calculations.basic_calculations import Operation

logger = logging.getLogger(__name__)
import sys

from collections import namedtuple

__author__ = 'zmbq'


class OurParser(ArgumentParser):
    def error(self, message):
        print("Error: %s" % message, file=sys.stderr)
        print("Enter csm --help for help", file=sys.stderr)
        sys.exit(2)


def _create_parser():
    parser = OurParser(usage="\ncsm type input_molecule output_file [additional arguments]", allow_abbrev=False)

    parser.add_argument('type',
                        # choices=c_symmetries + s_symmetries + ['cs', 'ci', 'ch'],
                        help='The type of operation: cs, ci, ch, cN, sN')
    parser.add_argument('input', help='Input molecule file')
    parser.add_argument('output', default='output.txt', help='Output file')

    # Optional arguments (their names start with --)

    # types of calculations (default is exact):
    calculation_type = parser.add_argument_group('Calculation Type (default is exact)')
    calculation_type.add_argument('--approx', action='store_true', default=False,
                                  help='use the approximate algorithm to estimate the CSM')
    calculation_type.add_argument('--trivial', action='store_true', default=False,
                                  help='CSM of identity perm, or, if chains, CSM of chain permutation with no atom permutation')


    parser.add_argument('--timeout', default=300,
                        help="Specify a timeout for CSM in seconds. Default is 5 minutes (300)", type=int)
    parser.add_argument('--sn-max', type=int, default=8, help='The maximal sn to try, relevant only for chirality')

    # general input/calculation arguments:
    # parser.add_argument('--ignore-hy', action='store_true', default=False, help='Ignore Hydrogen atoms in computations')
    input_type = parser.add_argument_group('Input Arguments')
    input_type.add_argument('--remove-hy', action='store_true', default=False,
                            help='Remove Hydrogen atoms, rebuild molecule without them, and compute')
    input_type.add_argument('--ignore-sym', action='store_true', default=False,
                            help='Ignore all atomic symbols, performing a purely geometric operation')
    input_type.add_argument('--use-mass', action='store_true', default=False,
                            help='Use the atomic masses to define center of mass')
    input_type.add_argument('--babel-bond', action='store_true', default=False, help='Let OpenBabel compute bonding')
    # parser.add_argument('--no-babel',  action='store_true', default=False, help='force suppress automatically using OpenBabel to compute bonds')
    input_type.add_argument('--use-sequence', action='store_true', default=False,
                            help='create equivalence class for pdb file using sequence information. Can\'t be used with --use-chains')
    input_type.add_argument('--ignore-chains', action='store_true', default=False,
                            help='When a molecule has chains, ignore them (affects trivial, approx)')
    input_type.add_argument('--read-fragments', action='store_true', default=False,
                            help='Read fragments from .mol or .pdb file as chains')


    # calculation arguments that only apply to exact:
    exact_args = parser.add_argument_group('Arguments for Exact Algorithm')
    exact_args.add_argument('--use-perm', type=str,
                            help='EXACT ONLY: Compute exact CSM for a single permutation')
    exact_args.add_argument('--keep-structure', action='store_true', default=False,
                            help='EXACT ONLY: Maintain molecule structure from being distorted in the exact calculation')
    exact_args.add_argument('--no-constraint', action='store_true', default=False,
                            help='EXACT ONLY: Do not use the constraints algorithm to traverse the permutation tree')

    # calculation arguments that only apply to approx
    # parser.add_argument('--use-dir', type=str,
    #                    help='Run the approx algorithm using predefined axes as the starting point')
    approx_args = parser.add_argument_group('Arguments for Approx Algorithm')
    approx_args.add_argument('--detect-outliers', action='store_true', default=False,
                             help="APPROX ONLY:Use outlier detection to improve guesses for initial directions in approx algorithm")
    approx_args.add_argument('--no-orthogonal', action='store_true', default=False,
                             help="APPROX ONLY:Don't add orthogonal directions to calculated directions")
    approx_args.add_argument('--use-best-dir', action='store_true', default=False,
                             help='APPROX ONLY:Only use the best direction')
    approx_args.add_argument('--many-chains', action='store_true', default=False,
                             help='APPROX ONLY: Use the new chains algorithm for many chains. Will automatically apply use-chains')
    approx_args.add_argument('--greedy', action='store_true', default=False,
                             help='APPROX ONLY: use the old greedy approx algorithm (no hungarian)')
    approx_args.add_argument('--dir', nargs=3, type=float,
                             help='run approximate algorithm using a specific starting direction')

    # output formatting and printing options
    out_args = parser.add_argument_group("Output Arguments")
    out_args.add_argument('--format', help='Use a specific input/output format')
    out_args.add_argument('--json-output', action='store_true', default=False,
                        help='Print output in json format to a file')
    out_args.add_argument('--print-local', action='store_true', default=False,
                        help='Print the local CSM (csm for each atom) in the output file')
    out_args.add_argument('--output-perms', action='store', default=None,
                        help='Writes all enumerated permutations to file')
    out_args.add_argument('--output-branches', action='store_true', default=False,
                        help='Writes all backtracking branches to the console')
    out_args.add_argument('--print-approx', action='store_true', default=False,
                        help='add some printouts to approx')

    # defunct: no longer applied in code
    # parser.add_argument('--no-limit', action='store_true', default=False, help='Allows running program while ignoring computational complexity')
    # parser.add_argument('--babel-test', action='store_true', default=False, help="Test if the molecule is legal or not")
    # parser.add_argument('--time-only', action='store_true', default=False, help="Only print the time and exit")
    # parser.add_argument('--write-openu', action='store_true', default=False,
    #                    help='Write output in open university format')


    return parser







def _process_arguments(parse_res):
    """
    Divides the parsed arguments (from argparse) into three dictionaries
    Args:
        parse_res: Result of argparse

    Returns:
        dictionary_args - input arguments, calculation arguments and output arguments

    """

    dictionary_args = {}

    # the first three positional arguments

    op = Operation(parse_res.type)

    dictionary_args['operation'] = op
    dictionary_args['op_type'] = op.type
    dictionary_args['op_order'] = op.order
    dictionary_args['op_name'] = op.name

    dictionary_args['in_file_name'] = parse_res.input

    dictionary_args['out_file_name'] = parse_res.output

    # optional arguments:
    dictionary_args['timeout'] = parse_res.timeout

    # types of calculations:
    dictionary_args['calc_type'] = 'exact'
    if parse_res.approx:
        dictionary_args['calc_type'] = 'approx'
    if parse_res.trivial:
        if parse_res.approx:
            raise ValueError("--approx and --trivial are mutually exclusive")
        if parse_res.just_perms:
            raise ValueError("--just-perms and --trivial are mutually exclusive")
        dictionary_args['calc_type'] = 'trivial'

    # general input/calculation arguments:
    # dictionary_args['ignore_hy'] = parse_res.ignore_hy
    dictionary_args['remove_hy'] = parse_res.remove_hy
    dictionary_args['ignore_symm'] = parse_res.ignore_sym
    dictionary_args['sn_max'] = parse_res.sn_max
    dictionary_args['use_mass'] = parse_res.use_mass
    dictionary_args['babel_bond'] = parse_res.babel_bond
    dictionary_args['use_sequence'] = parse_res.use_sequence

    # if parse_res.use_sequence and parse_res.keep_structure:
    #    raise ValueError("--keep-structure and --use-sequence are mutually exclusive")



    # use chains and fragments
    dictionary_args['use_chains'] = not parse_res.ignore_chains
    dictionary_args['read_fragments'] = parse_res.read_fragments

    if not dictionary_args['use_chains'] and parse_res.read_fragments:
        dictionary_args['use_chains'] = True
        logger.warn(
            "--read-fragments is only relevant when --use-chains has been specified, so --use-chains has been specified automatically")

    # calculation arguments for exact only:
    if parse_res.use_perm:
        if dictionary_args['calc_type'] != 'exact':
            logger.warning("--use-perm applies only to exact calculation.")
        dictionary_args['perm_file_name'] = parse_res.use_perm

    dictionary_args['keep_structure'] = parse_res.keep_structure
    if dictionary_args['calc_type'] in ['approx', 'trivial'] and parse_res.keep_structure:
        logger.warning("--keep-structure has no effect on approximate or trivial algorithms.")

    dictionary_args['no_constraint'] = parse_res.no_constraint
    if dictionary_args['calc_type'] in ['approx', 'trivial'] and parse_res.no_constraint:
        logger.warning("--no-constraint has no effect on approximate or trivial algorithms.")

    # calculation arguments for approx only:
    dictionary_args['approx_algorithm'] = 'hungarian'
    if parse_res.many_chains:
        if dictionary_args['calc_type'] != 'approx':
            logger.warning("--many-chains applies only to approx calculation. --many-chains will be ignored")
        if parse_res.greedy:
            raise ValueError("--many-chains and --greedy are mutually exclusive")
        dictionary_args['use_chains'] = True
        dictionary_args['approx_algorithm'] = 'many-chains'
    if parse_res.greedy:
        if dictionary_args['calc_type'] != 'approx':
            logger.warning("--greedy applies only to approx calculation. --greedy will be ignored")
        dictionary_args['approx_algorithm'] = 'greedy'

    dictionary_args['detect_outliers'] = parse_res.detect_outliers
    if dictionary_args['calc_type'] != 'approx' and parse_res.detect_outliers:
        logger.warning("--detect-outliers applies only to approx calculation. --detect-outliers will be ignored")

    dictionary_args['get_orthogonal'] = not parse_res.no_orthogonal
    if dictionary_args['calc_type'] != 'approx' and parse_res.no_orthogonal:
        logger.warning("--no-orthogonal applies only to approx calculation. --no-orthogonal will be ignored")

    dictionary_args['use_best_dir'] = parse_res.use_best_dir
    if dictionary_args['calc_type'] != 'approx' and parse_res.use_best_dir:
        logger.warning("--use-best-dir applies only to approx calculation. --use-best-dir will be ignored")

    dir = parse_res.dir
    if dir:
        dictionary_args['dirs'] = [dir]
    # if parse_res.use_dir:
    #    if dictionary_args['calc_type'] != 'approx':
    #        logger.warning("--use-dir applies only to approx calculation. --use-dir will be ignored")
    #    dictionary_args['dir_file_name'] = parse_res.use_dir




    # output arguments:
    dictionary_args['json_output'] = parse_res.json_output
    dictionary_args['print_approx'] = parse_res.print_approx
    dictionary_args['print_perms'] = parse_res.output_perms
    dictionary_args['print_branches'] = parse_res.output_branches
    dictionary_args['format'] = parse_res.format
    dictionary_args['useformat'] = dictionary_args['format'] is not None
    if not dictionary_args['format']:
        # get input file extension
        dictionary_args['format'] = parse_res.input.split(".")[-1]

    # dictionary_args['write_openu'] = parse_res.write_openu
    dictionary_args['print_local'] = dictionary_args['calc_local'] = parse_res.print_local

    dictionary_args['perms_csv_name'] = parse_res.output_perms

    permuters.print_branches = parse_res.output_branches

    return dictionary_args


def get_split_arguments(args):
    """
    :param args:
    :return:
    """
    parser = _create_parser()
    parsed_args, leftovers = parser.parse_known_args(args)
    dictionary_args = _process_arguments(parsed_args)
    return dictionary_args
