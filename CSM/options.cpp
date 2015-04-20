/*
* The CSM runtime options
*
* Extracted from mainRot.cpp during the 2014 reorganization
* Created by Itay Zandbank
*/

#include "options.h"
#include "logging.h"
#include <cstring>
#include <stdio.h>

using namespace std;

void csm_options::usage(const std::string op) 
{
	printf("Usage: %s <type> input_file output_file [-options]\n", op.c_str());
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
}

csm_options::csm_options() { }

csm_options::csm_options(int argc, char *argv[])
{
	char innerOpName[100];
	// check number of arguments
	if (argc < 4){
		usage(argv[0]);
		exit(1);
	}

	if (!strcmp(argv[1], "cs")) {
		type = CS;
		opOrder = 2;
		sprintf(innerOpName, "MIRROR SYMMETRY");
	}
	else if (!strcmp(argv[1], "ci")) {
		type = CI;
		opOrder = 2;
		sprintf(innerOpName, "INVERSION (S2)");
	}
	else if (!strcmp(argv[1], "ch")) {
		type = CH;
		opOrder = 2;
		sprintf(innerOpName, "CHIRALITY");
	}
	else if (argv[1][0] == 'c') {
		type = CN;
		opOrder = atoi(argv[1] + 1);
		sprintf(innerOpName, "C%d SYMMETRY", opOrder);
	}
	else if (argv[1][0] == 's') {
		type = SN;
		opOrder = atoi(argv[1] + 1);
		if (opOrder % 2 != 0) {
			LOG(fatal) << "Only Even values of n are allowed";
			exit(1);
		}
		sprintf(innerOpName, "S%d SYMMETRY", opOrder);
	}
	opName = innerOpName;

	// try to open infile for reading
	inFileName = argv[2];
	if ((inFile = fopen(inFileName.c_str(), "rt")) == NULL){
		if (writeOpenu)
		{
			printf("ERR* Failed to open data file %s *ERR\n", inFileName);
		}
		LOG(fatal) << "Failed to open data file " << inFileName;
		exit(1);
	}


	// try to open outfile for writing
	outFileName = argv[3];
	if ((outFile = fopen(outFileName.c_str(), "w")) == NULL){
		if (writeOpenu) {
			printf("ERR* Failed to open output file %s for writing *ERR\n", outFileName);
		}
		LOG(fatal) << "Failed to open output file " << outFileName << " for writing";
		exit(1);
	}

	// get commandline flags
	int i;
	int nextIsPermFile = false;
	int nextIsMaxSn = false;
	int nextIsDirFile = false;
	for (i = 4; i< argc; i++){
		if (nextIsPermFile) {
			char* permfileName = argv[i];
			if ((permfile = fopen(permfileName, "rt")) == NULL){
				if (writeOpenu)
				{
					printf("ERR* Failed to open perm file %s for reading *ERR\n", permfileName);
				}
				LOG(fatal) << "Failed to open perm file " << permfileName << "f or reading";
				exit(1);
			}
			nextIsPermFile = false;
		}
		else if (nextIsDirFile) {
			char* dirfilename = argv[i];
			if ((dirfile = fopen(dirfilename, "rt")) == NULL){
				if (writeOpenu) {
					printf("ERR* Failed to open dir file %s for reading *ERR\n", dirfilename);
				}
				LOG(fatal) << "Failed to open dir file " << dirfilename << " for reading";
				exit(1);
			}
			nextIsDirFile = false;
		}
		else if (nextIsMaxSn) {
			sn_max = atoi(argv[i]);
			nextIsMaxSn = false;
		}
		else if (strcmp(argv[i], "-sn_max") == 0) {
			if (type != CH) {
				LOG(fatal) << "Option -sn_max only applies to chirality";
				exit(1);
			}
			nextIsMaxSn = true;
		}
		else if (strcmp(argv[i], "-ignoreHy") == 0)
			ignoreHy = true;
		else if (strcmp(argv[i], "-removeHy") == 0)
			removeHy = true;

		else if (strcmp(argv[i], "-ignoreSym") == 0)
			ignoreSym = true;

		else if (strncmp(argv[i], "-format", 7) == 0) {
			useFormat = true;
			format = argv[i] + 7;
		}
		else if (strcmp(argv[i], "-writeOpenu") == 0) {
			writeOpenu = true;
		}
		else if (strcmp(argv[i], "-nolimit") == 0) {
			limitRun = false;
		}
		else if (strcmp(argv[i], "-useperm") == 0) {
			useperm = true;
			nextIsPermFile = true;
		}
		else if (strcmp(argv[i], "-usedir") == 0) {
			useDir = true;
			nextIsDirFile = true;
		}
		else if (strcmp(argv[i], "-babelbond") == 0) {
			babelBond = true;
		}
		else if (strcmp(argv[i], "-useMass") == 0) {
			useMass = true;
		}
		else if (strcmp(argv[i], "-timeonly") == 0) {
			timeOnly = true;
		}
		else if (strcmp(argv[i], "-printNorm") == 0) {
			printNorm = true;
		}
		else if (strcmp(argv[i], "-help") == 0) {
			usage(argv[0]);
			exit(0);
		}
		else if (strcmp(argv[i], "-findperm") == 0) {
			findPerm = true;
		}
		else if (strcmp(argv[i], "-detectOutliers") == 0) {
			detectOutliers = true;
		}
		else if (strcmp(argv[i], "-approx") == 0) {
			detectOutliers = true;
			findPerm = true;
		}
		else if (strcmp(argv[i], "-babelTest") == 0) {
			babelTest = true;
		}
		else if (strcmp(argv[i], "-printlocal") == 0) {
			printLocal = true;
		}
		else if (strcmp(argv[i], "-keepCenter") == 0) {
			keepCenter = true;
		}
		else if (strcmp(argv[i], "-log") == 0) {
			logFile = argv[i + 1];
			if (i + 1 == argc)
			{
				LOG(fatal) << "The -log option must be followed by a filename";
				exit(1);
			}
			i++;  // Skip the next argument
		}
	}
	if (writeOpenu) {
		useFormat = true;
		format = "PDB";
	}
}
