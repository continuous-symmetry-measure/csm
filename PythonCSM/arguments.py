"""
Parse the CSM command line arguments.

        Example:  CSM c2 molecule.xyz output.txt
        c2 is the type, molecule.xyz is the input_file, output.txt is the output file

        CSM s4 molecule.xyz output.txt -ignoreSym -useperm file.perm

        CSM s4     <--- Fails - not enough arguments

		printf("Usage: %s <type> input_file output_file [-options]\n",)
		printf("type is one of: cs, ci, cn, sn (n is replaced by the number between 2 and 9: c2. ... c9, s2, ... s9)\n");
		printf("Available options are:\n");
		printf("-ignoreHy - ignore Hydrogen atoms in computations\n");
		printf("-removeHy - remove Hydrogen atoms in computations,\n");
		printf("	rebuild molecule without them and compute\n");
		printf("-ignoreSym - Ignore all atomic symbols,\n");
		printf("	performing a purely geometric operation\n");
		printf("-formatXXX - Use a specific input/output format");
		printf("-writeOpenu - Write output in open university format\n");
		printf("-nolimit - Allows running program while ignoring computational complexity\n");
		printf("-useperm permfile  -  only compute for a single permutation\n");
		printf("-usedir dirfile  -  use a predefined axis as a starting point\n");
		printf("	This options ignores the -ignoreSym/-ignoreHy/-removeHy flags\n");
		printf("-findperm	   - Attempt to search for a permutation\n");
		printf("-detectOutliers - Use statistical methods to try and improve -findperm's results\n");
		//	printf("-anneal		   - Try to anneal the result\n");
		printf("-babelbond	   - Let openbabel compute bonding\n");
		printf("-useMass	   - Use the atomic masses to define center of mass\n");
		printf("-timeOnly	   - Only print the time and exit\n");
		printf("-babelTest	   - Test if the molecule is legal or not\n");
		printf("-sn_max	<max n> - The maximal sn to try, relevant only for chirality\n");
		printf("-printNorm		- Print the normalization factor as well\n");
		printf("-printlocal		- Print the local CSM (csm for each atom) in the output file\n");
		printf("-approx			- Equivalent to -detectOutliers -findperm together\n");
		printf("-keepCenter		- Do not change coordinates s.t. (0,0,0) corresponds to Center of Mass\n");
		printf("-log [logfile]	- Write a detailed log to logfile");
		printf("-help - print this help file\n");

"""
from argparse import ArgumentParser
from argparse import ArgumentTypeError
import re

__author__ = 'zmbq'


def create_parser():
    parser = ArgumentParser()

    # The first three positional arguments
    parser.add_argument('type', choices=('c2', 'c3', 'c4', 'c8', 's2', 's3', 's4', 's8', 'cs', 'ci'),
                        help='The type of operation')
    parser.add_argument('input', help='Input file')
    parser.add_argument('output', help='Output file')

    # Optional arguments (their names start with --)
    parser.add_argument('--ignoreHy', action='store_true', default=False, help='Ignore Hydrogen atoms in computations')
    parser.add_argument('--removeHy', action='store_true', default=False,
                        help='Remove Hydrogen atoms in computations, rebuild molecule without them and compute')
    parser.add_argument('--ignoreSym', action='store_true', default=False,
                        help='Ignore all atomic symbols, performing a purely geometric operation')

    # TODO: add possible format options
    # TODO: "--formatXXX" of "--format XXX" ???
    parser.add_argument('--format', choices=('XXX', 'ABC'), help='Use a specific input/output format')
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
    parser.add_argument('--anneal', action='store_true', default=False, help='Try to anneal the result')
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
                        help='Equivalent to -detectOutliers -findperm together')
    parser.add_argument('--keepCenter', action='store_true', default=False,
                        help='Do not change coordinates s.t. (0,0,0) corresponds to Center of Mass')
    parser.add_argument('--log', type=str, help='Write a detailed log to logfile')

    return parser

if __name__ == '__main__':
    # This code runs when you execute the script from PyCharm.
    parser = create_parser()

    args = ['--help']
    args = ['c2', 'abc.xyz', 'abc.output', '--useperm', 'permfile.dat']
    args = ['c4', 'abc.xyz', 'def.txt', '--ignoreHy', '--removeHy', '--format', 'ABC', '--sn_max', '5', '--log', 'log.txt']

    result = parser.parse_args(args)  # Call the parse

    print(result.type)
    print(result.input)
    print(result.output)
    print("Ignore Hydrogen: ", result.ignoreHy)
    print("Remove Hydrogen: ", result.removeHy)
    print("Ignore Symbols: ", result.ignoreSym)
    print("Format: ", result.format)
    print("Write Openu: ", result.writeOpenu)
    print("No limit: ", result.nolimit)
    if result.useperm:
        print("Use permutation: ", result.useperm)
    else:
        print("No permutation specified")
    if result.usedir:
        print("Use as a starting axis: ", result.usedir)
    else:
        print("No starting axis specified")
    print("Find permutation: ", result.findperm)
    print("Detect outliers: ", result.detectOutliers)
    print("Anneal the result: ", result.anneal)
    print("OpenBabel bonding: ", result.babelbond)
    print("Use atomic masses: ", result.useMass)
    print("Time only: ", result.timeOnly)
    print("Babel test: ", result.babelTest)
    if result.sn_max:
        print("Maximal sn: %d" % result.sn_max)
    else:
        print("No maximal sn specified")
    print("Print normalization factor: ", result.printNorm)
    print("Print local: ", result.printLocal)
    print("Approx: ", result.approx)
    print("Keep center: ", result.keepCenter)
    if result.log:
        print("Log: ", result.log)
    else:
        print("No log file specified")


    # Parsed arguments are in 'parser'
    print("\n\nHello from arguments.py")