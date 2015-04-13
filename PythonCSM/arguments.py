"""
Parse the CSM command line arguments.

		printf("Usage: %s <type> input_file output_file [-options]\n",)
		printf("type is one of: cs, ci, cn, sn (n is replaced by the number)\n");
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
__author__ = 'zmbq'

def create_parser():
    # Create the argparse parser
    pass

if __name__=='__main__':
    parser = create_parser()
    args = ['c2', 'abc.xyz', 'abc.output', '--useperm']


