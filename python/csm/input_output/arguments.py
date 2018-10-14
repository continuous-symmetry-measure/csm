"""
Parse the CSM command line arguments.
"""
from argparse import ArgumentParser
import logging

from csm.calculations.data_classes import Operation

logger = logging.getLogger(__name__)
import sys


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

    # Optional arguments

    # types of calculations (default is exact):
    calculation_type = parser.add_argument_group('Calculation Type (default is exact)')
    _calculation_type=calculation_type.add_mutually_exclusive_group()
    _calculation_type.add_argument('--approx', action='store_const', default='exact', const='approx', dest='calc_type',
                                  help='use the approximate algorithm to estimate the CSM')
    _calculation_type.add_argument('--trivial', action='store_const', const='trivial', dest='calc_type',
                                  help='CSM of identity perm, or, if chains, CSM of chain permutation with no atom permutation')


    parser.add_argument('--timeout', default=300,
                        help="Specify a timeout for CSM in seconds. Default is 5 minutes (300)", type=int)
    parser.add_argument('--sn-max', type=int, default=8, help='The maximal sn to try, relevant only for chirality')

    # general input/calculation arguments:
    input_type = parser.add_argument_group('Input Arguments')
    input_type.add_argument('--remove-hy', action='store_true', default=False,
                            help='Remove Hydrogen atoms, rebuild molecule without them, and compute')
    input_type.add_argument('--ignore-sym', action='store_true', default=False,
                            help='Ignore all atomic symbols, performing a purely geometric operation')
    input_type.add_argument('--use-mass', action='store_true', default=False,
                            help='Use the atomic masses to define center of mass')
    input_type.add_argument('--babel-bond', action='store_true', default=False, help='Let OpenBabel compute bonding')
    #input_type.add_argument('--use-sequence', action='store_true', default=False,
    #                        help='create equivalence class for pdb file using sequence information. Can\'t be used with --use-chains')
    input_type.add_argument('--use-chains', action='store_true', default=False,
                            help='When a molecule has chains, use them (affects trivial, approx)')

    # calculation arguments that only apply to approx
    approx_args = parser.add_argument_group('Arguments for Approx Algorithm')
    _approx_args=approx_args.add_mutually_exclusive_group()
    _approx_args.add_argument('--many-chains', action='store_const', default='hungarian', const='many-chains', dest='approx_algorithm',
                             help='APPROX ONLY: Use the new chains algorithm for many chains. Will automatically apply use-chains')
    _approx_args.add_argument('--greedy', action='store_const', const='greedy', dest='approx_algrithm',
                             help='APPROX ONLY: use the old greedy approx algorithm (no hungarian)')

    return parser

def _process_arguments(parse_res):
    """
    Divides the parsed arguments (from argparse) into three dictionaries
    Args:
        parse_res: Result of argparse

    Returns:
        dictionary_args - input arguments, calculation arguments and output arguments

    """

    dictionary_args = dict(vars(parse_res))

    # the first three positional arguments

    op = Operation(parse_res.type)

    dictionary_args['operation'] = op
    dictionary_args['op_type'] = op.type
    dictionary_args['op_order'] = op.order
    dictionary_args['op_name'] = op.name

    dictionary_args['in_file_name'] = parse_res.input
    dictionary_args['out_file_name'] = parse_res.output

    if dictionary_args["calc_type"]=="approx":
        dictionary_args["use_sequence"] = True

    elif dictionary_args["approx_algorithm"]!="hungarian":
        logger.warning("--many-chains and --greedy are only relevant with --approx and will be ignored")


    return dictionary_args


def get_split_arguments(args):
    """
    :param args:
    :return:
    """
    parser = _create_parser()
    parsed_args=parser.parse_args(args)
    #parsed_args, leftovers = parser.parse_known_args(args)
    dictionary_args = _process_arguments(parsed_args)
    return dictionary_args


if __name__ == '__main__':
    args=sys.argv[1:]
    get_split_arguments(args)
