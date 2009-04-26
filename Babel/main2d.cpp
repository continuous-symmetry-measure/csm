/*
* Author: shadi lahham, modified by Amir Zayit
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
#include "nrutil.h"
}

#include <openbabel/mol.h>
#include "babelAdapter.h"

#define CSMFORMAT "CSM"
#define MAXDOUBLE  100000000.0
#define MINDOUBLE  1e-8
#define GROUPSIZE_LIMIT 15
#define GROUPSIZE_FACTOR 1.32e11
#define APPROX_RUN_PER_SEC 8e4
#define TRUE 1
#define FALSE 0

typedef enum {
	CN,
	CS
} OperationType;

extern "C" {

// from tqli.c - nrbook
double tqli(double d[], double e[], int n, double **z);
// from tred2.c - nrbook
void tred2(double **a, int n, double d[], double e[]);
// from rpoly.c Jenkins-Traub Polynomial Solver
int rpoly(double *op, int degree, double *zeror, double *zeroi);

}

// function declarations
void printOutput(Molecule* m, double** outAtoms, double csm, double *dir, double dMin, FILE *out);
void printOutputPDB(Molecule* m, double** outAtoms, double csm, double *dir, double dMin, FILE *out);
void printOutputFormat(Molecule* m, OBMol& mol, double** outAtoms, double csm, double *dir, double dMin, FILE *out, char *fname);
void csmOperation(Molecule* m, double** outAtoms, int *optimalPerm, double* csm, double* dir, double* dMin, OperationType type);
void runSinglePerm(Molecule* m, double** outAtoms, int* perm, double* csm, double* dir, double* dMin, OperationType type);
void readPerm(FILE* permfile, int* perm, int size);
void initIndexArrays(Molecule* m, int* posToIdx, int* idxToPos);
double createSymmetricStructure(Molecule* m, double **outAtom, int *perm, double *dir, OperationType type, double dMin);

// global options
int ignoreHy = FALSE;
int removeHy = FALSE;
int ignoreSym = FALSE;
int useFormat = FALSE;
int writeOpenu = FALSE;
OperationType type;
OperationType chMinType;
int opOrder;
int useperm    = FALSE;
int useMass    = FALSE;
int limitRun = TRUE;
char *format = NULL;
int babelBond = FALSE;
int timeOnly = FALSE;

// file pointers
FILE* inFile = NULL;
FILE* outFile = NULL;
FILE* permfile = NULL;
char *inFileName = NULL;
char *outFileName = NULL;

char opName[100];

void usage(char *op) {
	printf("Usage: %s <type> input_file output_file [-options]\n", op);
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
	printf("	This options ignores the -ignoreSym/-ignoreHy/-removeHy flags\n");
	printf("-babelbond	   - Let openbabel compute bonding\n");
	printf("-useMass	   - Use the atomic masses to define center of mass\n");
	printf("-timeOnly	   - Only print the time and exit\n");	
	printf("-help - print this help file\n");
}

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
double numPermutations(Molecule *m, int operationOrder, OperationType t) {
	int i,j, k; 	
	double total = 1.0;
	for (i = 1; i <= m->_groupNum; i++) {
		int groupSize = getGroupSize(m, i);
		double temp = 0;
		double fact = factorial(groupSize);
		for (j = 0; j <= (groupSize / operationOrder); j++) {
			temp += fact/(factorial(groupSize - j * operationOrder )*pow((double)operationOrder, 	
				j)*factorial(j));
		}
		total *= (temp);
	}
	return total;
}

double totalNumPermutations(Molecule *m) {	
	return numPermutations(m, opOrder, type);	
}

char *getExtension(char *fname) {
	return strrchr(fname,'.') + 1;
}

/*
* parses the command line parameters
*/
void parseInput(int argc, char *argv[]){

	// check number of arguments
	if (argc < 4){
		usage(argv[0]);
		exit(1);
	}

	if (!strcmp(argv[1],"cs")) {
		type = CS;
		opOrder = 2;
		sprintf(opName, "MIRROR SYMMETRY");			
	} else if (argv[1][0] == 'c') {
		type = CN;
		opOrder = atoi(argv[1] + 1);
		sprintf(opName, "C%d SYMMETRY",opOrder);	
	} 
	// try to open infile for reading
	inFileName = argv[2];
	if ((inFile = fopen(inFileName, "rt")) == NULL){
		if (writeOpenu) {
			printf("ERR* Failed to open data file %s *ERR\n", inFileName);
		} else {
			printf("Failed to open data file %s\n", inFileName);
		}
		exit(1);
	}


	// try to open outfile for writing
	outFileName = argv[3];
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
	int nextIsMaxSn = FALSE;

	for ( i=4;  i< argc ;  i++ ){
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
		} else if (strcmp(argv[i],"-ignoreHy" ) == 0 )
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
		} else if (strcmp(argv[i], "-useMass") == 0) { 
			useMass = TRUE;					
		} else if (strcmp(argv[i], "-timeonly") == 0) {
			timeOnly = TRUE;
		} else if (strcmp(argv[i], "-help") == 0) {
			usage(argv[0]);
			exit(0);
		}
	}
	if (writeOpenu) {
		useFormat = TRUE;
		format = strdup("PDB");		
	}
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
	int *perm = NULL;	

	// init options
	parseInput(argc,argv);

	// try to read molecule from infile
	Molecule* m;
	OBMol mol; 
	if (useFormat) {
		// If a specific format is used, read molecule using that format
		if (strcasecmp(format, CSMFORMAT) == 0) {
			m = createMolecule(inFile,stdout,ignoreSym && !useperm);
			if (useMass) { 
				for (int i = 0; i < m->_size; i++) { 
					m->_mass[i] = getAtomicMass(m->_symbol[i]);
				}
			} 
		} else {
			mol = readMolecule (inFileName, format, babelBond);
			m = babel2Mol(mol, ignoreSym && !useperm, useMass);	
		}
   	} else {
		format = getExtension(inFileName);

		// if the extension is CSM - use csm
		if (strcasecmp(format, CSMFORMAT) == 0) {
			m = createMolecule(inFile,stdout,ignoreSym && !useperm);
			if (useMass) { 
				for (int i = 0; i < m->_size; i++) { 
					m->_mass[i] = getAtomicMass(m->_symbol[i]);
				}
			}
		} else {
			
			mol = readMolecule (inFileName, NULL, babelBond);
			m = babel2Mol(mol, ignoreSym && !useperm, useMass);			
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

	printf("Going to enumerate over %5.2f permutations\n", totalNumPermutations(m));
	printf("Entire run should take approx. %5.2f hours on a 2.0Ghz Computer\n", 1.0*totalNumPermutations(m) / 3600 / 		APPROX_RUN_PER_SEC );
	if (timeOnly) { return 0; };

	// allocate memory for outAtoms
	outAtoms = (double **)malloc(m->_size * sizeof(double*));
	for (i=0;i<m->_size;i++)
		outAtoms[i] = (double *)malloc(3 * sizeof(double));
       
	perm = (int *)malloc(sizeof(int) * m->_size);

	//normalize Molecule
	if (!normalizeMolecule(m)){
		if (writeOpenu) {
			printf("ERR* Failed to normalize atom positions: dimension of set of points = zero *ERR\n");
		} else {
			printf("Failed to normalize atom positions: dimension of set of points = zero\n");
		}
		exit(1);
	}

	if (useperm) {	
		readPerm(permfile,perm, m->_size);
		runSinglePerm(m, outAtoms, perm, &csm, dir, &dMin, type);
	} else {
		// perform operation
		csmOperation(m, outAtoms, perm, &csm, dir, &dMin, type);			
	}
	
	// Denormalize
	for (i = 0; i < m->_size; i++) { 
		m->_pos[i][0] *= m->_norm;
		m->_pos[i][1] *= m->_norm;
		m->_pos[i][2] *= m->_norm;
		outAtoms[i][0] *= m->_norm;
		outAtoms[i][1] *= m->_norm;
		outAtoms[i][2] *= m->_norm;
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

	fprintf(outFile, "\n PERMUTATION:\n\n");
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

/**
 * Calculate the best axis, and compute the CSM for it, given a pairing of the indices (perm)
 */ 
double calcRefPlane(Molecule* m, int* perm, double *dir, OperationType type) {

	if (type == CN) { 
		int *curPerm = (int *)malloc(sizeof(int) * m->_size);
		double rotated[2], dists, csm;
		int i, j;
	
		// initialize identity permutation
		for (i = 0; i < m->_size; i++) {
			curPerm[i] = i;
		}
	
		// compute CSM	
		csm = 0.0;		
		for (i = 0; i < opOrder; i++) {		
			if (i != 0) {
				double angle  = (2 * PI * i) / opOrder;
				dists = 0.0;
				double c = cos(angle);
				double s = sin(angle);
				for (j = 0; j < m->_size; j++) {
					// i'th power of permutation
					curPerm[j] = perm[curPerm[j]];	
				}
				
				for (j = 0; j < m->_size; j++) {
					rotated[0] = c * m->_pos[curPerm[j]][0] - s * m->_pos[curPerm[j]][1];
					rotated[1] = s * m->_pos[curPerm[j]][0] + c* m->_pos[curPerm[j]][1];
					dists += (m->_pos[j][0] * rotated[0] + m->_pos[j][1] * rotated[1]); 
				}
				csm += dists;
			} else {
				csm += 1.0;	
			}
		}	
		// DEBUG printf("%.6f %.6f %.6f\n", csm, maxval, scl);
		// DEBUG printf("dir: %.6f %.6f %.6f\n", dir[0], dir[1], dir[2]);
		// DEBUG printf("eig: %.6f %.6f %.6f\n",  diag[1], diag[2], diag[3]);
		csm = fabs(100 * (1.0 - csm / opOrder));
		free(curPerm);

	} else if (type == CS) {
	}
	return csm;
}

double createSymmetricStructure(Molecule* m, double **outAtoms, int *perm, double *dir, OperationType type, double dMin) {


	if (type == CN) { 
		int *curPerm = (int *)malloc(sizeof(int) * m->_size);
		double rotated[2];
		int i,j,k;
		double res;
	
		for (i = 0; i < m->_size; i++) {
			// initialize with identity operation
			curPerm[i] = i;
			for (j = 0; j < 2; j++) {
				outAtoms[i][j] = m->_pos[i][j];
			}
		}
	
		// compute CSM	
		for (i = 0; i < opOrder; i++) {		
			if (i != 0) {
				double angle  = (2 * PI * i) / opOrder;
				double c = cos(angle);
				double s = sin(angle);
				for (j = 0; j < m->_size; j++) {
					// i'th power of permutation
					curPerm[j] = perm[curPerm[j]];	
				}
				
				for (j = 0; j < m->_size; j++) {
					outAtoms[j][0] += c * m->_pos[curPerm[j]][0] - s * m->_pos[curPerm[j]][1];
					outAtoms[j][1] += s * m->_pos[curPerm[j]][0] + c* m->_pos[curPerm[j]][1];
				}
			}
		}	
		free(curPerm);


	
		for (j = 0; j < m->_size; j++) {
			for (k = 0; k < 2; k++) {
				outAtoms[j][k] /= opOrder;
				outAtoms[j][k] *= dMin;
				res += SQR(outAtoms[j][k]);
			}
		}
		return sqrt(res);
	} else if (type == CS) {
		return -1;
	} 
}


/*
* Calculates minimal csm, dMin and directional cosines by applying the chiralityFunction
* breaks down Molecule into groups of similar atoms and calculates the above only for
* permutations that keep the similar atoms within the group ( groupPermuter class )
* once it finds the optimal permutation , calls the chiralityFunction on the optimal permutation
*/
void csmOperation(Molecule* m, double** outAtoms, int *optimalPerm, double* csm, double* dir, 
			double* dMin, OperationType type){

	int i;
	double curCsm;
	double curDir[3];
	int *idxToPos, *posToIdx;
	int * groupSizes;
	groupPermuter* gp;
	int addGroupsOfTwo;

	// These are the 
	int *realPerm = (int *)malloc(sizeof(int) * m->_size);

	// allocate memory for index arrays arrays
	idxToPos = (int*)malloc(m->_size * sizeof(int));
	posToIdx = (int*)malloc(m->_size * sizeof(int));

	// and group sizes arrays
	groupSizes = (int*)malloc(m->_groupNum * sizeof(int*));	

	for(i=0; i<m->_size; i++)    	
		// init index arrays
		initIndexArrays(m,posToIdx,idxToPos);

	// init csm, curCsm
	*csm = curCsm = MAXDOUBLE;

	// get group sizes
	for(i=1; i<= m->_groupNum ; i++){
		groupSizes[i-1] = getGroupSize(m,i);
	}

	
	addGroupsOfTwo = 0;	
	gp = createGroupPermuter(m->_groupNum,groupSizes,m->_size,opOrder, addGroupsOfTwo);
	if (!gp){
		if (writeOpenu) {
			printf("ERR* Failed to create groupPermuter *ERR\n");
		} else {
			printf("Failed to create groupPermuter \n");	
		}	
		exit(1);
	};

	// calculate csm for each valid permutation & remember minimal (in optimalAntimer)
	while ( nextGroupPermutation(gp) ) {

		for (i = 0; i < m->_size; i++) {
			realPerm[i] = idxToPos[gp->_index[posToIdx[i]]];			
		}			
		curCsm = calcRefPlane(m, realPerm, curDir, type);		

		// check, if it's a minimal csm, update maxGroupCsm and optimalAntimer
		if(curCsm < *csm) {
			*csm = curCsm;
			dir[0] = curDir[0]; dir[1] = curDir[1]; dir[2] = curDir[2];
			
			for (i = 0; i < m->_size; i++) {
				optimalPerm[i] = realPerm[i];
			}
		}				
	}

	// failed to find value for any permutation
	
	if (*csm == MAXDOUBLE){
		if (writeOpenu) {
			printf("ERR* Failed to calculate a csm value for %s *ERR\n",opName);
		} else {
			printf("Failed to calculate a csm value for %s \n",opName);
		}
		exit(1);
	}

	// which is DMIN?
	*dMin = (1.0 - (*csm / 100 * opOrder / (opOrder - 1)));
	createSymmetricStructure(m, outAtoms, optimalPerm, dir, type, *dMin);
		

	// housekeeping
	free(groupSizes);
	free(idxToPos);
	free(posToIdx);
	freeGroupPermuter(gp);
	free(realPerm);
}

/*
 * Calculates csm, dMin and directional cosines for a given permutation
 */
void runSinglePerm(Molecule* m, double** outAtoms, int *perm, double* csm, double* dir, double* dMin, OperationType type){
	*csm = calcRefPlane(m, perm, dir, type);


	// which is DMIN?
	*dMin = (1.0 - (*csm / 100 * opOrder / (opOrder - 1)));
	createSymmetricStructure(m, outAtoms, perm, dir, type, *dMin);
}



/*
* prints the Molecule position, outcome position, csm, dMin and directional cosines to output file
*/
void printOutput(Molecule* m, double** outAtoms, double csm, double *dir, double dMin, FILE *out){

	int i,j;
	printf("%s: %.6lf\n",opName,fabs(csm));
	fprintf(out, "%s: %.4lf\n",opName,fabs(csm));
	fprintf(out, "SCALING FACTOR: %7lf\n", dMin);

	fprintf(out, "\n INITIAL STRUCTURE COORDINATES %i\n\n",m->_size);

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

	fprintf(out, "\n RESULTING STRUCTURE COORDINATES %i\n\n",m->_size);

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

	fprintf(out, "\n DIRECTIONAL COSINES:\n\n");
	fprintf(out, "%lf %lf %lf\n", dir[0], dir[1], dir[2]);

}

/*
* prints PDB ATOM tags
*/
void printPDBATOM(Molecule* m,FILE* f,char** sym,double** pos){
	int i;
	for(i=0; i<m->_size; i++){
		fprintf(f,"ATOM  %5d %2s                %8.3lf%8.3lf%8.3lf                      %2s\n",
			i+1,sym[i],pos[i][0],pos[i][1],pos[i][2],sym[i]);
	}
}

/*
* prints PDB CONECT tags
*/
void printPDBCONNECT(Molecule* m,FILE* f){
	int i,j;
	for(i=0; i<m->_size; i++){
		fprintf(f,"CONECT%5d",i +1);
		for ( j=0;  j< m->_valency[i] ; j++ ){
			if ((j>0) && (!(j%4)))
				fprintf(f,"\nCONECT%5d",i +1);
			fprintf(f,"%5d",m->_adjacent[i][j] +1);
		}
		fprintf(f,"\n");
	}
}

/*
* prints in PDB format the Molecule position, outcome position, csm, dMin and directional cosines to output file
*/
void printOutputPDB(Molecule* m, double** outAtoms, double csm, double *dir, double dMin, FILE *out){

	//	out = stdout;
	// print PDB file
	fprintf(out,"MODEL        1\n");

	printPDBATOM(m,out,m->_symbol,m->_pos);

	printPDBCONNECT(m,out);

	fprintf(out,"ENDMDL\n");
	fprintf(out,"MODEL        2\n");

	printPDBATOM(m,out,m->_symbol,outAtoms);

	printPDBCONNECT(m,out);

	fprintf(out,"ENDMDL\n");

	// print results to screen


	if (writeOpenu)
		printf("SV* %.4lf *SV\n",fabs(csm));
	else
		printf( "%s: %.4lf\n",opName,fabs(csm));
	printf( "SCALING FACTOR: %7lf\n", dMin);
	printf( "DIRECTIONAL COSINES: %lf %lf %lf\n", dir[0], dir[1], dir[2]);

}

/*
* prints in PDB format the Molecule position, outcome position, csm, dMin and directional cosines to output file
*/
void printOutputFormat(Molecule* m, OBMol& mol, double** outAtoms, double csm, double *dir, double dMin, FILE *out, char *fname) {

	fprintf(out, "%s: %.4lf\n",opName,fabs(csm));
	fprintf(out, "SCALING FACTOR: %7lf\n", dMin);

	// TODO - should we print the centered molecule, or the original one (and, accordingly, the symmetric struct)

	fprintf(out, "\n INITIAL STRUCTURE COORDINATES %i\n\n",m->_size);

	updateCoordinates(mol, m->_pos);

	writeMolecule(mol, format, out, fname);

	updateCoordinates(mol, outAtoms);

	fprintf(out, "\n RESULTING STRUCTURE COORDINATES %i\n\n",m->_size);

	writeMolecule(mol, format, out, fname);

	// print results to screen

	fprintf(out, "\n DIRECTIONAL COSINES:\n\n");
	fprintf(out, "%lf %lf %lf\n", dir[0], dir[1], dir[2]);

	if (writeOpenu)
		printf("SV* %.4lf *SV\n",fabs(csm));
	else
		printf( "%s: %.4lf\n",opName,fabs(csm));
	printf( "SCALING FACTOR: %7lf\n", dMin);
	printf( "DIRECTIONAL COSINES: %lf %lf %lf\n", dir[0], dir[1], dir[2]);

}




