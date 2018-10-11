from argparse import ArgumentParser, SUPPRESS
import logging

import os

from datetime import datetime

from csm.calculations import permuters
from csm.calculations.data_classes import Operation
from csm import __version__
logger = logging.getLogger(__name__)
import sys

class OurParser(ArgumentParser):
    def error(self, message):
        sys.stdout.write("Error: %s" % message)
        sys.stdout.write("\nEnter csm --help for help, or csm [command] --help for help with a specific command")
        sys.exit(2)

def _create_parser():

    def input_utility_func(parser):
        parser.add_argument('--connect', const=os.path.join(os.getcwd(), "connectivity.txt"), nargs='?',
                            help='xyz connectivity file, default is connectivity.txt in current working directory')
        parser.add_argument('--remove-hy', action='store_true', default=False,
                            help='Remove Hydrogen atoms, rebuild molecule without them, and compute')
        parser.add_argument('--ignore-sym', action='store_true', default=False,
                            help='Ignore all atomic symbols, performing a purely geometric operation')
        parser.add_argument('--use-mass', action='store_true', default=False,
                            help='Use the atomic masses to define center of mass')
        parser.add_argument('--babel-bond', action='store_true', default=False,
                            help='Let OpenBabel compute bonding')
        parser.add_argument('--use-sequence', action='store_true', default=False,
                            help='create equivalence class for pdb file using sequence information.')
        parser.add_argument('--use-chains', action='store_true', default=False,
                            help='When a molecule has chains, use them (affects trivial, approx)')
        parser.add_argument('--read-fragments', action='store_true', default=False,
                            help='Read fragments from .mol or .pdb file as chains')
        parser.add_argument('--select-mols', default=None,
                            help='Select only some molecules, eg 1-20,15,17,19-21')
        parser.add_argument('--select-atoms',default=None,
                            help='Select only some atoms, eg 1-20,15,17,19-21')

    def output_utility_func(parser):
        parser.add_argument('--json-output', action='store_true', default=False,
                            help='Print output in json format to a file. Only relevant with --legacy')
        parser.add_argument('--print-local', action='store_true', default=False,
                            help='Print the local CSM (csm for each atom) in the output file')
        parser.add_argument('--print-denorm', action='store_true', default=False,
                            help='when printing the original molecule, print the denormalized coordinates')
        parser.add_argument("--legacy",  action='store_true', default=False,
                            help='print the old csm format')
        parser.add_argument("--simple",  action='store_true', default=False,
                            help='only output is CSM to screen')
        parser.add_argument("--not-unique", action='store_true', default=False,
                            help="don't ensure results folder is unique")

    def shared_calc_utility_func(parser):
        parser.add_argument('type',
                            # choices=c_symmetries + s_symmetries + ['cs', 'ci', 'ch'],
                            help='The type of operation: cs, ci, ch, cN, sN',)
                            #nargs="+")
        parser.add_argument('--timeout', default=300,
                            help="Specify a timeout for CSM in seconds. Default is 5 minutes (300)", type=int)
        parser.add_argument('--global-timeout', default=None,
                            help="Specify a global timeout for CSM in seconds. Default is 50000 seconds (over 13 hours)", type=int)
        parser.add_argument('--sn-max', type=int, default=8,
                            help='The maximal sn to try, relevant only for chirality')
        parser.add_argument("--pipe", action='store_true', default=False,
                            help="treat this program as a piped program (read from sys.stdin, write to sys.stdout)")

    def shared_normalization_utility_func(parser): #I made this because having normalization stuck in the calc utility func was ugly
        parser.add_argument('--normalize', default=[],
                            help='Types of normalization available:\n'
                                 '0: standard normalization, according to centers of mass (without scaling)\n'
                                 '1: normalization according to the center of mass of each fragment\n'
                                 '2: normalization according to an approximation of the symmetric structure of the centers '
                                 'of mass of each fragment, based on the solution permutation\n'
                                 '3: normalization according to an approximation of the symmetric structure of the centers '
                                 'of mass of each fragment, without using the solution permutation\n'
                                 '4: normalization according to averages of approximation to symmetry of fragments\n'
                                 '5: normalization according to number of atoms\n'
                                 '6: linear normalization',
                            choices=['0', '1', '2', '3', '4', '5', '6'],
                            nargs='+'
                            )

    def add_input_output_utility_func(parser):
        parser_input_args = parser.add_argument_group("Args for input (requires --input)")
        parser_input_args.add_argument("--input", help="molecule file or folder, default is current working directory", const=os.getcwd(), nargs='?')
        parser_input_args.add_argument('--in-format', help='override guessing format from input file ending with provided format',
                                default=None)
        input_utility_func(parser_input_args)
        parser_output_args = parser.add_argument_group("Args for output (requires --output)")
        parser_output_args.add_argument("--output", const=os.path.join(os.getcwd(), 'csm_results', timestamp), nargs='?',
                                        help="output file or folder, default is 'csm_results\\timestamp' folder in current working directory, if provided directory exists a new one with timestamp will be created",)
        parser_output_args.add_argument('--out-format', help='override guessing format from output file ending with provided format',
                                default=None)
        output_utility_func(parser_output_args)


    parser = OurParser(allow_abbrev=False, usage="csm read/write/exact/trivial/approx/comfile [args] \n"
                                                 "example: csm exact c2 --input mymol.mol --output --keep-structure\n"
                                                 "for specific help with each subprogram and its available arguments, enter csm COMMAND -h\n"
                                                 "e.g. csm exact -h")
    timestamp = str(datetime.now().timestamp())[-11:].replace(".", "")
    parser.add_argument('--timestamp', help=SUPPRESS, default=timestamp)
    parser.add_argument("--version", help="print version and exit", action='store_true', default=False)
    commands = parser.add_subparsers(title="Available commands", dest="command")

    #command
    commands_args_=commands.add_parser('comfile', help='provide a command file for running calculations')
    command_args=commands_args_.add_argument_group("Command args")
    command_args.add_argument('comfile', default=os.path.join(os.getcwd(), "cmd.txt"), nargs='?',
                              help="the file that contains the commands, default is cmd.txt in current working directory")
    command_args.add_argument('--old-cmd', action='store_true', default=False,
                              help="the old format with csm sym __INPUT__ __OUTPUT__ --approx etc")
    command_args.add_argument("--verbose", action='store_true', default=False,
                             help='create a fixed width spreadsheet of information about each direction for ann approx commands')
    add_input_output_utility_func(commands_args_)

    #READ
    input_args = commands.add_parser('read', help="Read a molecule file into a json in CSM format",
                                     usage="csm read filename [optional args]\n"
                                           "example: csm read mymol.pdb --read-fragments --remove-hy --select-atoms 1-3")
    input_args.add_argument('input', help='molecule file or folder, default is current working directory', default=os.getcwd(), nargs='?')
    input_args.add_argument('--format', help='override guessing format from file ending with provided format', default=None)
    input_utility_func(input_args)

    #WRITE
    out_args = commands.add_parser('write', help="Output the results of the calculation to a file- must be used with piped input",
                                   usage="csm write filename [optional args]")
    out_args.add_argument('output', default=os.path.join(os.getcwd(), 'csm_results', timestamp), nargs='?',
                          help="output file or folder, default is 'csm_results\\timestamp' folder in current working directory, if provided directory exists a new one with timestamp will be created",)
    out_args.add_argument('--format', help='override guessing format from file ending with provided format',
                            default=None)
    output_utility_func(out_args)

    #EXACT
    exact_args_ = commands.add_parser('exact', help="Perform an exact CSM calculation", conflict_handler='resolve',
                                      usage='csm exact TYPE [optional args]\n'
                                            'example: csm exact s4 --input --output myresults/1 --keep-structure --timeout 500')
    exact_args = exact_args_.add_argument_group("Args for exact calculation")
    shared_calc_utility_func(exact_args)
    exact_args.add_argument('--use-perm', nargs="?", type=str, default=None, const=os.path.join(os.getcwd(), "perm.txt"),
                            help='Compute exact CSM for a single permutation, default is current directory/perm.txt')
    exact_args.add_argument('--keep-structure', action='store_true', default=False,
                            help="Don't allow permutations that break bonds")
    exact_args.add_argument('--output-branches', action='store_true', default=False,
                            help="Don't allow permutations that break bonds")
    exact_args.add_argument('--output-perms', const="DEFAULT", nargs='?',
                            help='Writes all enumerated permutations to file. Default is OUTPUT_DIR/perms.csv, or working directory/perms.csv is --output not selected')
    shared_normalization_utility_func(exact_args)
    add_input_output_utility_func(exact_args_)


    #APPROX
    approx_args_ = commands.add_parser('approx', help="Approximate the CSM value", conflict_handler='resolve',
                                       usage='csm approx TYPE [optional args]\n'
                                             'example: csm approx ch --input --output --detect-outliers --parallel 4 --sn-max 10')
    approx_args = approx_args_.add_argument_group("Args for approx calculation")
    shared_calc_utility_func(approx_args)
    #choosing dir:
    approx_args.add_argument('--detect-outliers', action='store_true', default=False,
                             help="Use outlier detection to improve guesses for initial directions in approx algorithm")
    approx_args.add_argument('--no-orthogonal', action='store_true', default=False,
                             help="Don't add orthogonal directions to calculated directions")
    approx_args.add_argument('--fibonacci', type=int,
                             help="Use fibonacci sphere to generate N starting directions")
    approx_args.add_argument('--use-best-dir', action='store_true', default=False,
                             help='Only use the best direction')
    approx_args.add_argument('--dir', nargs=3, type=float,
                             help='run approximate algorithm using a specific starting direction')
    #algorithm choice
    approx_args.add_argument('--greedy', action='store_true', default=False,
                             help='use the old greedy approx algorithm (no hungarian)')
    approx_args.add_argument('--many-chains', action='store_true', default=False,
                             help='Use the new chains algorithm for many chains. Will automatically apply use-chains')
    approx_args.add_argument('--keep-structure', action='store_true', default=False,
                             help='Use keep-structure approximate algorithm')
    approx_args.add_argument('--selective', type=int,
                             help='Do a single iteration on many directions (use with --fibonacci), and then a full set of iterations only on the best k (default 10)')
    approx_args.add_argument('--parallel', type=int, const=0, nargs='?',
                             help='Calculate directions in parallel. Recommended for use with fibonacci. If no number of processors is specified, cpu count - 1 will be used')
    #outputs
    approx_args.add_argument("--verbose", action='store_true', default=False,
                             help='create a fixed width spreadsheet of information about each direction')
    approx_args.add_argument('--polar', action='store_true', default=False,
                             help="Print polar coordinates instead of cartesian coordinates in statistics")
    approx_args.add_argument('--print-approx', action='store_true', default=False,
                             help='print log to screen from approx')
    shared_normalization_utility_func(approx_args)
    add_input_output_utility_func(approx_args_)


    #TRIVIAL
    trivial_args_ = commands.add_parser('trivial', help="Calculate trivial (identity) CSM", conflict_handler='resolve',
                                        usage='csm trivial TYPE [optional args]\n'
                                              'example: csm trivial c4 --input --output --permute-chains')
    trivial_args = trivial_args_.add_argument_group("Args for trivial calculation")
    shared_calc_utility_func(trivial_args)
    #this is totally equivalent to --use-chains, however --use-chains is under input arguments and I want permute chains to have
    #documentation specifically under calculation arguments for trivial, as it's THE main calculation choice for trivial
    trivial_args.add_argument('--permute-chains', action='store_true', default=False, help="Permute the chains before calculating trivial calculation. "
                                                                                           "Will automatically activate --use-chains, and is automatically activated if --use-chains is provided")
    shared_normalization_utility_func(trivial_args)
    add_input_output_utility_func(trivial_args_)
    return parser

def _process_arguments(parse_res):
    def _parse_ranges_and_numbers(input):
        '''
        given a string such as '1-5,8,12-13', where 1 is the first index, returns an array of indices that matches
        (but where 0 is the first index) (1-5 includes 5)
        :param input: a string with number ranges and individual numbers, separated by commas
        :return: [int, int, int..]
        '''
        selected=[]
        if input:
            items=input.split(',')
            for item in items:
                if "-" in item:
                    range_limits = item.split("-")
                    for i in range(int(range_limits[0]), int(range_limits[1])+1):
                        selected.append(i - 1)
                else:
                    selected.append(int(item) - 1)
        return selected

    def parse_input(dictionary_args):
        dictionary_args['in_file_name'] = parse_res.input
        dictionary_args["conn_file"]=parse_res.connect
        if parse_res.read_fragments:
            dictionary_args['use_chains'] = True
            #logger.warning("Warning: --read-fragments is only relevant when --use-chains has been specified, so --use-chains has been specified automatically")

        dictionary_args['select_mols'] = _parse_ranges_and_numbers(parse_res.select_mols)
        dictionary_args['select_atoms']  = _parse_ranges_and_numbers(parse_res.select_atoms)
        if parse_res.select_atoms and parse_res.remove_hy:
            raise ValueError("Remove-hy and select-atoms cannot be used together")


    def parse_output(dictionary_args):
        dictionary_args['out_file_name'] = parse_res.output
        dictionary_args['print_local'] = dictionary_args['calc_local'] = parse_res.print_local

    dictionary_args=vars(parse_res)
    if "pipe" not in dictionary_args:
        dictionary_args["pipe"]=False

    if parse_res.command == "read":
        dictionary_args["in_format"]=parse_res.format
        parse_input(dictionary_args)
    elif parse_res.command == "write":
        dictionary_args["out_format"] = parse_res.format
        parse_output(dictionary_args)
    else:
        # get input/output if relevant
        parse_input(dictionary_args)
        parse_output(dictionary_args)

        if parse_res.command == "comfile":
            dictionary_args["command_file"]=parse_res.comfile
            dictionary_args["old_command"]=parse_res.old_cmd
        else:
            dictionary_args['operation'] = Operation(parse_res.type)
            dictionary_args['normalizations']=parse_res.normalize

            if parse_res.command == 'exact':
                if parse_res.use_perm:
                    dictionary_args['perm_file_name'] = parse_res.use_perm
                dictionary_args['print_branches'] = parse_res.output_branches
                permuters.print_branches = parse_res.output_branches

                dictionary_args['perms_csv_name'] = parse_res.output_perms
                if parse_res.output_perms=="DEFAULT":
                    try:
                        base_path= os.path.dirname(os.path.abspath(dictionary_args["out_file_name"]))
                        dictionary_args['perms_csv_name']=os.path.join(base_path,  "perms.csv")
                    except:
                        dictionary_args['perms_csv_name'] = os.path.join("perms.csv")
            if parse_res.command == 'approx':
                #choose dir:
                #dictionary_args['detect_outliers'] = parse_res.detect_outliers
                dictionary_args['get_orthogonal'] = not parse_res.no_orthogonal
                if parse_res.fibonacci is not None:
                    dictionary_args["num_dirs"] = parse_res.fibonacci
                    dictionary_args["fibonacci"] = True #doing this before previous line causes weird bug
                dir = parse_res.dir
                if dir:
                    dictionary_args['dirs'] = [dir]

                #algorithm choice:
                dictionary_args['approx_algorithm'] = 'hungarian'
                if parse_res.greedy:
                    dictionary_args['approx_algorithm'] = 'greedy'
                if parse_res.many_chains:
                    if parse_res.greedy:
                        raise ValueError("--many-chains and --greedy are mutually exclusive")
                    dictionary_args['use_chains'] = True
                    dictionary_args['approx_algorithm'] = 'many-chains'
                if parse_res.keep_structure:
                    dictionary_args['approx_algorithm'] = 'structured'

                if parse_res.selective is not None:
                    dictionary_args["num_selected"] = parse_res.selective
                    dictionary_args["selective"] = True #likewise, weird bugs


                if parse_res.parallel is not None:
                    dictionary_args['pool_size'] = parse_res.parallel
                    dictionary_args['parallel'] = True #likewise, weird bugs

                    #outputs:
                #dictionary_args['print_approx'] = parse_res.print_approx
                #dictionary_args['polar'] = parse_res.polar
            if parse_res.command == 'trivial':
                if parse_res.permute_chains or parse_res.use_chains:
                    dictionary_args["use_chains"]=True

    try:
        timestamp=str(parse_res.timestamp)
        out_file_name=dictionary_args['out_file_name']
        if os.path.exists(out_file_name) and not parse_res.not_unique:
                if not os.path.isfile(out_file_name):
                    head, tail = os.path.split(out_file_name)
                    dictionary_args['out_file_name'] = os.path.join(head, tail + timestamp)
                else:
                    filename=os.path.basename(out_file_name)
                    filename=filename[:-4] + timestamp + filename[-4:]
                    head, tail = os.path.split(out_file_name)
                    dictionary_args['out_file_name']=os.path.join(head, filename)

    except (KeyError, TypeError, AttributeError) as e:
    #there is no output, eg in Read
    #output is None, eg in a line of command
    #parse_res doesn't have unique attribute
        pass

    return dictionary_args

def get_parsed_args(args):
    parser = _create_parser()
    parsed_args = parser.parse_args(args)
    if parsed_args.version:
        print("CSM version:", __version__)
        sys.exit()
    if parsed_args.command is None:
        parser.error("You must select a command from: read, exact, approx, trivial, write")
    processed_args = _process_arguments(parsed_args)
    return processed_args



def get_allowed_args_for_command(command):
    parser = _create_parser()
    commands = parser._subparsers._group_actions[0]._name_parser_map
    allowed_args = commands[command]._option_string_actions
    return allowed_args

def check_modifies_molecule(cmd):
    modifies_molecule=False
    allowed_mol_args=get_allowed_args_for_command('read')
    for arg in cmd.split():
        if arg in allowed_mol_args:
            modifies_molecule = True
            break
    return modifies_molecule

def old_cmd_converter(cmd):
    '''
    receives a command string. returns a valid (in the new args format) set of args.
    used by --command to read lines from a file
    :param cmd:
    :return:
    '''

    if cmd[:3]=="csm":
        cmd=cmd[3:]
    args = cmd.split()
    symm = args[0]
    input = args[1]
    output = args[2]
    command = "exact"
    if "--trivial" in args:
        command = "trivial"
    if "--approx" in args:
        command = "approx"

    final_args = [command, symm]
    allowed_args=get_allowed_args_for_command(command)

    prev_arg_fine = False
    for arg in args[3:]:
        if arg[:2] != "--":
            if prev_arg_fine:
                final_args.append(arg)
        if arg in allowed_args:
            final_args.append(arg)
            prev_arg_fine = True


        else:
            prev_arg_fine = False

    return final_args