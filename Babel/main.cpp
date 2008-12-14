
/*
 * Author: shadi lahham
 *
 * Main body that initiates chirality operations.
 *
 * deals with input output, and the main logic of the high level calculation
 *
 */

extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> //for strcmp,strlen
#include "Molecule.h"
#include "groupPermuter.h"
#include "mainhelpers.h"
}

#include <openbabel/mol.h>
#include "babelAdapter.h"

#define CSMFORMAT "CSM"
#define MAXDOUBLE  100000000.0
#define MINDOUBLE  1e-8
#define GROUPSIZE_LIMIT 15
#define GROUPSIZE_FACTOR 1.32e11
#define APPROX_RUN_PER_SEC 2e5
#define TRUE 1
#define FALSE 0

// function declarations
void printOutput(Molecule* m, double** outAtoms, double csm, double *dir, double dMin, FILE *out);
void printOutputPDB(Molecule* m, double** outAtoms, double csm, double *dir, double dMin, FILE *out);
void chiralityOperation(Molecule* m, double** outAtoms, double* csm, double* dir, double* dMin, int* bestPerm);
void runSinglePerm(Molecule* m, double** outAtoms, double* csm, double* dir, double* dMin, int* perm);
void readPerm(FILE* permfile, int* perm, int size);
void initIndexArrays(Molecule* m, int* posToIdx, int* idxToPos);
void printOutputFormat(Molecule* m, OBMol& mol, double** outAtoms, double csm, double *dir, double dMin, FILE *out, char *fname);

// global options
int ignoreHy = FALSE;
int removeHy = FALSE;
int ignoreSym = FALSE;
int useFormat = FALSE;
int writeOpenu = FALSE;
int useperm    = FALSE;
int limitRun = TRUE;
char *format = NULL;
int babelBond = FALSE;
int timeOnly = FALSE;

// file pointers
FILE* inFile = NULL;
FILE* outFile = NULL;
FILE* permfile = NULL;
char* inFileName;
char* outFileName;

void usage(char *op) {
	printf("Usage: %s input_file output_file [-options]\n", op);
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
	printf("	This options ignores the -ignoreSym/-ignoreHy/-removeHy flags\n");
	printf("-babelbond	   - Let openbabel compute bonding\n");
	printf("-timeonly	   - Quit after printing time\n");
	printf("-help - print this help file\n");
}

char *getExtension(char *fname) {
	return strrchr(fname,'.') + 1;
}

/*
 * parses the command line parameters
 */
void parseInput(int argc, char *argv[]){

	// check number of arguments
	if (argc < 3){
		usage(argv[0]);
		exit(1);
	}

	// try to open infile for reading
    inFileName = argv[1];
    if ((inFile = fopen(inFileName, "rt")) == NULL){
		if (writeOpenu) {
			printf("ERR* Failed to open data file %s *ERR\n", inFileName);
		} else {
			printf("Failed to open data file %s\n", inFileName);
		}
        exit(1);
    }

    // try to open outfile for writing
    outFileName = argv[2];
    if ((outFile = fopen(outFileName, "w")) == NULL){
		if (writeOpenu) {
	        printf("ERR* Failed to open output file %s for writing *ERR\n", outFileName);
		} else {
			printf("Failed to open output file %s for writing\n", outFileName);
		}
        exit(1);
    }

    // get commandline flags
	int i;
	int nextIsPermFile = FALSE;

	for ( i=3;  i< argc ;  i++ ){
		if (nextIsPermFile) {
			char* permfileName = argv[i];
			if ((permfile = fopen(permfileName, "rt")) == NULL){
				if (writeOpenu) {
					printf("ERR* Failed to open perm file %s for reading *ERR\n", permfileName);
				} else {
					printf("Failed to open perm file %s for reading\n", permfileName);
				}
				exit(1);
			}
			nextIsPermFile = FALSE;
		}
	    else if (strcmp(argv[i],"-ignoreHy" ) == 0 )
	    	ignoreHy = TRUE;

	    else if (strcmp(argv[i],"-removeHy" ) == 0 )
	    	removeHy = TRUE;

	    else if (strcmp(argv[i],"-ignoreSym" ) == 0 )
	    	ignoreSym = TRUE;

	    else if (strncmp(argv[i],"-format", 7 ) == 0 ) {
		useFormat = TRUE;
		format = strdup(argv[i] + 7);
	    } else if (strcmp(argv[i],"-writeOpenu" ) == 0 ) {
	    	writeOpenu = TRUE; 	    
	    } else if (strcmp(argv[i], "-nolimit") == 0) {
			limitRun = FALSE; 
		} else if (strcmp(argv[i], "-useperm") == 0) {
			useperm = TRUE;
			nextIsPermFile = TRUE;
		} else if (strcmp(argv[i], "-babelbond") == 0) {
			babelBond = TRUE;			
		} else if (strcmp(argv[i], "-help") == 0) {
			usage(argv[0]);
			exit(0);
		} else if (strcmp(argv[i], "-timeonly") == 0) {
			timeOnly = TRUE;
		} 
	}
	if (writeOpenu) {
		useFormat = TRUE;
		format = strdup("PDB");		
	}
}

// ************************************************************
//       helper functions
// ************************************************************

/*
 * reutnrs n!
 */
double factorial(int n){
	double fact = 1.0;
  	while(n > 0)
    	fact *= n--;
    return fact;
}

/*
 * checks if the molecule's group permutations exceed GROUPSIZE_FACTOR_PER_MINUTE * MINUTES
 */
double numPermutations(Molecule *m){
	int i,j; 
	double total = 1.0;
	for (i = 1; i <= m->_groupNum; i++) {
		int groupSize = getGroupSize(m, i);
		double temp = 0;
		double fact = factorial(groupSize);
		for (j = 0; j <= (groupSize / 2); j++) {
			temp += fact/(factorial(groupSize - j * 2 )*pow(2.0, j)*factorial(j));
		}
		total *= (temp);
	}
	return total;
}

// ************************************************************
//       main program
// ************************************************************

/*
 * main funciton - check valid parameters, parse molecule and call chirality Operation
 */
int main(int argc, char *argv[]){

   int i;
   double csm, dMin;
   double **outAtoms;                 // output atom coordinates
   double dir[3] = {0.0, 0.0, 0.0};   // directional cosines
   int *perm = NULL;				   // optimal (or single) permutation

   // init options
   parseInput(argc,argv);

   // try to read molecule from infile
   Molecule* m;
   OBMol mol;
   if (useFormat) {
	// If a specific format is used, read molecule using that format
	if (strcasecmp(format, CSMFORMAT) == 0) {
		m = createMolecule(inFile,stdout,ignoreSym && !useperm);
	} else {
		mol = readMolecule (inFileName, format, babelBond);
		m = babel2Mol(mol, ignoreSym && !useperm );	
	}
   } else {
	format = getExtension(inFileName);	

	// if the extension is CSM - use csm
	if (strcasecmp(format, CSMFORMAT) == 0) {
		m = createMolecule(inFile,stdout,ignoreSym && !useperm);
	} else {
		mol = readMolecule (inFileName, NULL, babelBond);
		m = babel2Mol(mol, ignoreSym && !useperm );	
	}
   }

   if (!m){
		if (writeOpenu) {
	        printf("ERR* Failed to read molecule from data file *ERR\n");
		} else {
			printf("Failed to read molecule from data file \n");
		}
        exit(1);
    }

    // strip unwanted atoms if needbe
    if ((ignoreHy || removeHy) && !useperm){
		char* removeList[] = {"H"," H"};
		Molecule* n = NULL;
		if (ignoreHy)
			n = stripAtoms(m,removeList,2,FALSE);
		else //removeHy
			n = stripAtoms(m,removeList,2,TRUE);		
		mol.DeleteHydrogens();		

		if (!n){
			if (writeOpenu) {
	        	printf("ERR* Failed while trying to strip unwanted atoms *ERR\n");
			} else {
				printf("Failed while trying to strip unwanted atoms \n");
			}
        	exit(1);
    	}
		freeMolecule(m);
		m = n;
		// continue as per usual
    }

	if (!writeOpenu) {
		printf("Going to enumerate over %5.2f permutations\n", numPermutations(m));
		printf("Entire run should take approx. %5.2f hours on a 2.0Ghz Computer\n", 1.0*numPermutations(m) / 3600 / 	APPROX_RUN_PER_SEC );
		if (timeOnly) { return 0; };
	}

	if (limitRun && !useperm) {
		// check if groupsize exceeds limits
		if (getMaxGroupSize(m) > GROUPSIZE_LIMIT){
			if (writeOpenu) {
		        printf("ERR* Molecule contains subgroups of a size exceeding the limit of %d *ERR\n", GROUPSIZE_LIMIT);
			} else {
				printf("Molecule contains subgroups of a size exceeding the limit of %d\n", GROUPSIZE_LIMIT);
			}
			exit(1);
		}

		// check if the permutation of groupsizes is too complex
		if (numPermutations(m) > GROUPSIZE_FACTOR){
			if (writeOpenu) {
		        printf("ERR* Molecule structure is too complex and would require a long computation time *ERR\n");
			} else {
				printf("Molecule structure is too complex and would require a long computation time \n");
			}
			exit(1);
		}
	}

    	// allocate memory for outAtoms
	outAtoms = (double **)malloc(m->_size * sizeof(double*));
	for (i=0;i<m->_size;i++)
		outAtoms[i] = (double *)malloc(3 * sizeof(double));

	// perform operation
	perm = (int *)malloc(sizeof(int) * m->_size);
	if (useperm) {
		readPerm(permfile,perm, m->_size);
		runSinglePerm(m, outAtoms, &csm, dir, &dMin,perm);
	} else {
		chiralityOperation(m, outAtoms, &csm, dir, &dMin, perm);
	}

   	if (useFormat) {
		// If a specific format is used, read molecule using that format
		if (strcasecmp(format, CSMFORMAT) == 0) {
			printOutput(m, outAtoms, csm, dir, dMin, outFile);
		} else {
			printOutputFormat(m, mol, outAtoms, csm, dir, dMin, outFile, outFileName);
		}
	} else {
		// if the extension is CSM - use csm
		if (strcasecmp(getExtension(inFileName), CSMFORMAT) == 0) {
			printOutput(m, outAtoms, csm, dir, dMin, outFile);
		} else {
			printOutputFormat(m, mol, outAtoms, csm, dir, dMin, outFile, outFileName);
		}
	}

	fprintf(outFile, "\n PERMUTATIION:\n\n");
	for (i = 0; i < m->_size; i++) {
		fprintf(outFile, "%d ", perm[i] + 1);
	}
	fprintf(outFile,"\n");

	// housekeeping
	for (i=0;i<m->_size;i++){
		free(outAtoms[i]);
	}
	free(outAtoms);
    	freeMolecule(m);
	free(perm);

    	fclose(inFile);
    	fclose(outFile);
	if (permfile != NULL)
		fclose(permfile);
    	return 0;
}

/**
 * Read a permutation from a file - simply a list separated by spaces.
 */
void readPerm(FILE* permfile, int* perm, int size) {
	int *used = (int *)malloc(sizeof(int) * size);
	int i = 0;
	for (i = 0; i < size; i++) {
		used[i] = FALSE;
	}
	for (i = 0; i < size; i++) {
		int cur = -1;
		int res = fscanf(permfile,"%d", &cur);
		if (res != 1 || cur < 1 || cur > size || used[cur - 1]) {
			if (writeOpenu) {
				printf("ERR* Invalid permutation *ERR\n");
			} else {
				printf("Invalid permutation\n");
			}
			free(used);
			fclose(permfile);
			exit(1);
		}
		used[cur - 1] = TRUE;
		perm[i] = cur - 1;
	}
	free(used);
}


/*
 * creates position to index , and index to position translation arrays
 */
void initIndexArrays(Molecule* m, int* posToIdx, int* idxToPos){

	int i,j,counter;

	counter = 0;

	// build idxToPos
	for ( j=1;  j<= m->_groupNum ;  j++ ){
		for ( i=0;  i< m->_size ;  i++ ){
			if (m->_similar[i] == j){
				idxToPos[counter] = i;
				counter++;
			}
		}
	}

	// build posToIdx
	for ( i=0;  i< m->_size ;  i++ ){
		posToIdx[idxToPos[i]] = i;
	}

}

/*
 * Calculates minimal csm, dMin and directional cosines by applying the chiralityFunction
 * breaks down Molecule into groups of similar atoms and calculates the above only for
 * permutations that keep the similar atoms within the group ( groupPermuter class )
 * once it finds the optimal permutation , calls the chiralityFunction on the optimal permutation
 */
void chiralityOperation(Molecule* m, double** outAtoms, double* csm, double* dir, double* dMin,int *bestPerm){

	int i;
	double curCsm;
	int *idxToPos, *posToIdx;
	int * groupSizes;
	// antimer is a permutation of the position of atoms (or their mirror)
	double** antimer, **curAntimer, **optimalAntimer;
	double eulerParam[4];	

	//normalize Molecule
    if (!normalizeMolecule(m)){
    	printf("Failed to normalize atom positions: dimension of set of points = zero\n");
		exit(1);
    }

	// allocate memory for index arrays arrays
    idxToPos = (int*)malloc(m->_size * sizeof(int));
    posToIdx = (int*)malloc(m->_size * sizeof(int));

    // and for antimers
    curAntimer = (double**)malloc(m->_size * sizeof(double*));
    optimalAntimer = (double**)malloc(m->_size * sizeof(double*));

	// and group sizes arrays
	groupSizes = (int*)malloc(m->_groupNum * sizeof(int*));

    // init optimalAntimer - will get updated further down but just in case
    for(i=0; i<m->_size; i++)
    	optimalAntimer[i] = NULL;

	// init index arrays
	initIndexArrays(m,posToIdx,idxToPos);

    // init csm, curCsm
    *csm = curCsm = MAXDOUBLE;

    // create antimer
    antimer = createAntimer(m);
    

    // get group sizes
	for(i=1; i<= m->_groupNum ; i++){
		groupSizes[i-1] = getGroupSize(m,i);
	}
	for (i = 0; i < m->_size; i++) {
		bestPerm[i] = -1;
	}

	// create a groupPermuter
	extern int operationOrder;
	groupPermuter* gp = createGroupPermuter(m->_groupNum,groupSizes,m->_size,operationOrder, 0);
    if (!gp){
		if (writeOpenu) {
	        printf("ERR* Failed to create groupPermuter *ERR\n");
		} else {
			printf("Failed to create groupPermuter \n");	
		}
        exit(1);
    }

    // calculate csm for each valid permutation & remember minimal (in optimalAntimer)
    while ( nextGroupPermutation(gp) ){

    	// permute curAntimer
    	for ( i=0;  i< m->_size ;  i++ ){
			curAntimer[i] = antimer[idxToPos[gp->_index[posToIdx[i]]]];
		}

		// calculate current csm
		curCsm = chiralityFunction(m->_pos, curAntimer, m->_size, dir, dMin, eulerParam);		

		// check, if it's a minimal csm, update maxGroupCsm and optimalAntimer
    		if(curCsm < *csm){
			
    			*csm = curCsm;
			for ( i=0;  i< m->_size ;  i++ ) {
				optimalAntimer[i] = curAntimer[i];
				bestPerm[i] = idxToPos[gp->_index[posToIdx[i]]];
			}
		}

		// exit loop if reached zero
		if (curCsm < MINDOUBLE) {		
			break;
		}

	}

	// failed to find value for any permutation
	extern char* operationName;
	if (*csm == MAXDOUBLE){
		if (writeOpenu) {
	        printf("ERR* Failed to calculate a csm value for %s *ERR\n",operationName);
		} else {
			printf("Failed to calculate a csm value for %s \n",operationName);
		}
        exit(1);
	}

	// get csm for the optimalAntimer - this is the csm we're looking for
	*csm = chiralityFunction(m->_pos, optimalAntimer, m->_size, dir, dMin, eulerParam);

	// correct csm -  Chaim Dryzun measure correction
	// old // *csm = 100.0/2.0*(1.0-sqrt(1.0-*csm));
	// formula: S = 100 * (n-1) / n * (1-sqrt(1-csm))
	//*csm = 100.0 * (operationOrder - 1) / operationOrder * (1.0-sqrt(1.0-*csm));
	*csm = 100 * (operationOrder - 1) / operationOrder * (*csm);	

	
	foldingUnfolding(m->_pos, optimalAntimer, outAtoms, m->_size, dir, *dMin, eulerParam);

	// housekeeping
	free(groupSizes);
	free(idxToPos);
	free(posToIdx);
	free(curAntimer);
	free(optimalAntimer);

	//causes a crash on linux, ok on windows. Must investigate. disabled for now
	//freeGroupPermuter(gp);

	freeAntimer(m,antimer);

}

/*
* Calculates csm, dMin and directional cosines for a given permutation
*/
void runSinglePerm(Molecule* m, double** outAtoms, double* csm, double* dir, double* dMin, int *perm){

	int i;

	// antimer is a permutation of the position of atoms (or their mirror)
	double** antimer, **optimalAntimer;
	double eulerParam[4];

	extern int operationOrder;

	//normalize Molecule
	if (!normalizeMolecule(m)){
		if (writeOpenu) {
			printf("ERR* Failed to normalize atom positions: dimension of set of points = zero *ERR\n");
		} else {
			printf("Failed to normalize atom positions: dimension of set of points = zero\n");
		}
		exit(1);
	}

	// and for antimers
	optimalAntimer = (double**)malloc(m->_size * sizeof(double*));

	// init optimalAntimer - will get updated further down but just in case
	for(i=0; i<m->_size; i++) {
		optimalAntimer[i] = NULL;
	}

	// init csm, curCsm
	*csm = MAXDOUBLE;

	// create antimer
	antimer = createAntimer(m);

	// calculate csm for each valid permutation & remember minimal (in optimalAntimer)
	// permute curAntimer
	for ( i=0;  i< m->_size ;  i++ ){
		optimalAntimer[i] = antimer[perm[i]];
	}

	// get csm for the optimalAntimer - this is the csm we're looking for
	*csm = chiralityFunction(m->_pos, optimalAntimer, m->_size, dir, dMin, eulerParam);

	// correct csm -  Chaim Dryzun measure correction
	// old // *csm = 100.0/2.0*(1.0-sqrt(1.0-*csm));
	// formula: S = 100 * (n-1) / n * (1-sqrt(1-csm))
	//*csm = 100.0 * (operationOrder - 1) / operationOrder * (1.0-sqrt(1.0-*csm));
	*csm = 100 * (operationOrder - 1) / operationOrder * (*csm);

	foldingUnfolding(m->_pos, optimalAntimer, outAtoms, m->_size, dir, *dMin, eulerParam);

	// housekeeping
	free(optimalAntimer);

	//causes a crash on linux, ok on windows. Must investigate. disabled for now
	//freeGroupPermuter(gp);

	freeAntimer(m,antimer);

}


/*
 * prints the Molecule position, outcome position, csm, dMin and directional cosines to output file
 */
void printOutput(Molecule* m, double** outAtoms, double csm, double *dir, double dMin, FILE *out){

	int i,j;
	extern char* operationName;
	printf("%s: %.4lf\n",operationName,fabs(csm));
	fprintf(out, "%s: %.4lf\n",operationName,fabs(csm));
	fprintf(out, "SCALING FACTOR: %7lf\n", dMin);

	fprintf(out, "\n INITIAL STRUCTURE COORDINATES \n\n%i\n",m->_size);

	for(i=0; i<m->_size; i++){
		fprintf(out, "%3s%10lf %10lf %10lf\n",
		m->_symbol[i], m->_pos[i][0], m->_pos[i][1], m->_pos[i][2]);
	}

	for (i = 0; i < m->_size; i++) {
		fprintf(out, "%d ", i + 1);
		for ( j = 0; j < m->_valency[i]; j++ ) {
			fprintf(out, "%d ", m->_adjacent[i][j] + 1);
		}
		fprintf(out,"\n");
	}

	fprintf(out, "\n RESULTING STRUCTURE COORDINATES \n\n%i\n",m->_size);

	for(i=0; i<m->_size; i++){
		fprintf(out, "%3s%10lf %10lf %10lf\n",
		m->_symbol[i], outAtoms[i][0], outAtoms[i][1], outAtoms[i][2]);
	}

	for (i = 0; i < m->_size; i++) {
		fprintf(out, "%d ", i + 1);
		for ( j = 0; j < m->_valency[i]; j++ ) {
			fprintf(out, "%d ", m->_adjacent[i][j] + 1);
		}
		fprintf(out,"\n");
	}

	fprintf(out, "\n DIRECTIONAL COSINUSES:\n\n");
	fprintf(out, "%lf %lf %lf\n", dir[0], dir[1], dir[2]);

}

/*
* prints in PDB format the Molecule position, outcome position, csm, dMin and directional cosines to output file
*/
void printOutputFormat(Molecule* m, OBMol& mol, double** outAtoms, double csm, double *dir, double dMin, FILE *out, char *fname) {
	extern char* operationName;
	fprintf(out, "%s: %.4lf\n",operationName,fabs(csm));
	fprintf(out, "SCALING FACTOR: %7lf\n", dMin);

	fprintf(out, "\n INITIAL STRUCTURE COORDINATES %i\n\n",m->_size);

	writeMolecule(mol, format, out, fname);

	updateCoordinates(mol, outAtoms);

	fprintf(out, "\n RESULTING STRUCTURE COORDINATES %i\n\n",m->_size);

	writeMolecule(mol, format, out, fname);

	// print results to screen

	fprintf(out, "\n DIRECTIONAL COSINUSES:\n\n");
	fprintf(out, "%lf %lf %lf\n", dir[0], dir[1], dir[2]);

	if (writeOpenu)
		printf("SV* %.4lf *SV\n",fabs(csm));
	else
		printf( "%s: %.4lf\n",operationName,fabs(csm));
	printf( "SCALING FACTOR: %7lf\n", dMin);
	printf( "DIRECTIONAL COSINUSES: %lf %lf %lf\n", dir[0], dir[1], dir[2]);

}




