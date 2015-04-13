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

__author__ = 'zmbq'

def create_parser():
    parser = ArgumentParser()

    # The first three positional arguments
    parser.add_argument('type', choices=('c2', 'c3', 'c4', 'c8', 's2', 's3', 's4', 's8', 'cs', 'ci'),
                        help='The type of operation')
    parser.add_argument('input', help='Input file')
    # Add output

    # Optional arguments (their names start with --
    parser.add_argument('--ignoreHy', action='store_true', default=False, help='Ignore Hydrogen atoms in computations')
    parser.add_argument('--sn_max', type=int, help='The maximal sn to try, relevant only for chirality')
    return parser

if __name__=='__main__':
    # This code runs when you execute the script from PyCharm.
    parser = create_parser()

    # args = ['--help']
    args = ['c2', 'abc.xyz', 'abc.output', '--useperm', 'permfile.dat']
    args=['c4', 'abc.xyz', '--ignoreHy', '--sn_max', '5']
    result = parser.parse_args()  # Call the parse

    print(result.type)
    print(result.input)
    print("Ignore Hydrogen: ", result.ignoreHy)

    if result.sn_max:
        print("sn_max: %d" % result.sn_max)
    else:
        print("No sn_max specified")

    # Parsed arguments are in 'parser'
    print("Hello from arguments.py")


