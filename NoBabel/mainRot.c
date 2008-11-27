/*
* Author: shadi lahham, modified by Amir Zayit
*
* Main body that initiates chirality operations.
*
* deals with input output, and the main logic of the high level calculation
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> //for strcmp,strlen
#include "Molecule.h"
#include "groupPermuter.h"
#include "mainhelpers.h"
#include "nrutil.h"

#define MAXDOUBLE  100000000.0
#define MINDOUBLE  1e-8
#define GROUPSIZE_LIMIT 15
#define GROUPSIZE_FACTOR 1.32e11
#define APPROX_RUN_PER_SEC 2e5
#define TRUE 1
#define FALSE 0
#define SN_MAX 4

typedef enum {
	CN,
	SN, 
	CS,
	CI, 
	CH,
} OperationType;

// from tqli.c - nrbook
double tqli(double d[], double e[], int n, double **z);
// from tred2.c - nrbook
void tred2(double **a, int n, double d[], double e[]);
// from rpoly.c Jenkins-Traub Polynomial Solver
int rpoly(double *op, int degree, double *zeror, double *zeroi);

// function declarations
void printOutput(Molecule* m, double** outAtoms, double csm, double *dir, double dMin, FILE *out);
void printOutputPDB(Molecule* m, double** outAtoms, double csm, double *dir, double dMin, FILE *out);
void csmOperation(Molecule* m, double** outAtoms, int *optimalPerm, double* csm, double* dir, double* dMin, OperationType type);
void runSinglePerm(Molecule* m, double** outAtoms, int* perm, double* csm, double* dir, double* dMin, OperationType type);
void readPerm(FILE* permfile, int* perm, int size);
void initIndexArrays(Molecule* m, int* posToIdx, int* idxToPos);
double createSymmetricStructure(Molecule* m, double **outAtom, int *perm, double *dir, OperationType type, double dMin);

// global options
int ignoreHy = FALSE;
int removeHy = FALSE;
int ignoreSym = FALSE;
int readPDB = FALSE;
int writePDB = FALSE;
int writeOpenu = FALSE;
int useperm    = FALSE;
OperationType type;
OperationType chMinType;
int opOrder;
int timeOnly = FALSE;

// file pointers
FILE* inFile = NULL;
FILE* outFile = NULL;
FILE* permfile = NULL;

char opName[100];

void usage(char *op) {
	printf("Usage: %s <type> input_file output_file [-options]\n", op);
	printf("type is one of: cs, ci, cn, sn (n is replaced by the number\n");
	printf("Available options are:\n");
	printf("-ignoreHy - ignore Hydrogen atoms in computations\n");
	printf("-removeHy - remove Hydrogen atoms in computations,\n");
	printf("	rebuild molecule without them and compute\n");
	printf("-ignoreSym - Ignore all atomic symbols,\n");
	printf("	performing a purely geometric operation\n");
	printf("-readPDB - Use PDB as the input format, instead of XYZ\n");
	printf("-writePDB - Use PDB as the output format, instead of XYZ\n");
	printf("-writeOpenu - Write output in open university format\n");
	printf("-nolimit - Allows running program while ignoring computational complexity\n");
	printf("-useperm permfile  -  only compute for a single permutation\n");
	printf("	This options ignores the -ignoreSym/-ignoreHy/-removeHy flags\n");
	printf("-timeOnly	   - Only print the time and exit\n");
	printf("-help - print this help file");
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
double numPermutations(Molecule *m, int operationOrder) {
	int i,j; 
	double total = 1.0;
	for (i = 1; i <= m->_groupNum; i++) {
		int groupSize = getGroupSize(m, i);
		double temp = 0;
		double fact = factorial(groupSize);
		for (j = 0; j <= (groupSize / operationOrder); j++) {
			temp += fact/(factorial(groupSize - j * operationOrder )*pow(operationOrder, j)*factorial(j));
		}
		total *= (temp);
	}
	return total;
}

double totalNumPermutations(Molecule *m) {
	if (type != CH) { 
		return numPermutations(m, opOrder);	
	} else {
		// CS
		int i;
		double numPerms = numPermutations(m, 2);
		for (i = 2; i <= SN_MAX; i+=2) {
			numPerms += numPermutations(m, i);
		}	
		return numPerms;		
	}
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
	} else if (!strcmp(argv[1], "ci")) {
		type = CI;
		opOrder = 2;
		sprintf(opName, "INVERSION (S2)");	
	} else if (!strcmp(argv[1], "ch")) {
		type = CH;
		opOrder = 2;
		sprintf(opName, "CHIRALITY");	
	} else if (argv[1][0] == 'c') {
		type = CN;
		opOrder = atoi(argv[1] + 1);
		sprintf(opName, "C%d SYMMETRY",opOrder);	
	} else if (argv[1][0] == 's') {
		type = SN;
		opOrder = atoi(argv[1] + 1);
		if (opOrder % 2 != 0) {
			printf("ERROR - Only Even values of n are allowed\n");
			exit(1);
		}
		sprintf(opName, "S%d SYMMETRY",opOrder);	
	}
	
	// try to open infile for reading
	char* inFileName = argv[2];
	if ((inFile = fopen(inFileName, "rt")) == NULL){
		if (writeOpenu) {
			printf("ERR* Failed to open data file %s *ERR\n", inFileName);
		} else {
			printf("Failed to open data file %s\n", inFileName);
		}
		exit(1);
	}


	// try to open outfile for writing
	char* outFileName = argv[3];
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
		}
	    else if (strcmp(argv[i],"-ignoreHy" ) == 0 )
			ignoreHy = TRUE;

		else if (strcmp(argv[i],"-removeHy" ) == 0 )
			removeHy = TRUE;

		else if (strcmp(argv[i],"-ignoreSym" ) == 0 )
			ignoreSym = TRUE;

		else if (strcmp(argv[i],"-readPDB" ) == 0 )
			readPDB = TRUE;

		else if (strcmp(argv[i],"-writePDB" ) == 0 )
			writePDB = TRUE;

		else if (strcmp(argv[i],"-writeOpenu" ) == 0 )
			writeOpenu = TRUE; 
		else if (strcmp(argv[i], "-useperm") == 0) {
			useperm = TRUE;
			nextIsPermFile = TRUE;
		} else if (strcmp(argv[i], "-timeonly") == 0) {
			timeOnly = TRUE;
		} else if (strcmp(argv[i], "-help") == 0) {
			usage(argv[0]);
			exit(0);
		}
	}
	if (writeOpenu) {
		readPDB = TRUE;
		writePDB = TRUE;
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
	if (readPDB)
    	m = createMoleculePDB(inFile,stdout,ignoreSym && !useperm);
	else
    	m = createMolecule(inFile,stdout,ignoreSym && !useperm);
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

	// perform operation
	perm = (int *)malloc(sizeof(int) * m->_size);
	if (useperm) {	
		if (type == CH) {
			printf("Chirality can't be given a permutation, run the specific csm operation instead\n");
			exit(1);
		}	
		readPerm(permfile,perm, m->_size);
		runSinglePerm(m, outAtoms, perm, &csm, dir, &dMin, type);
	} else {
		if (type != CH) { 
			// perform operation
			csmOperation(m, outAtoms, perm, &csm, dir, &dMin, type);			
		} else {
			// chirality support 
			double chCsm, chdMin;
			double **chOutAtoms;
			int *chPerm;			
			double chDir[3] = {0.0, 0.0, 0.0};						

			chOutAtoms = (double **)malloc(m->_size * sizeof(double*));
			for (i=0;i<m->_size;i++)
				chOutAtoms[i] = (double *)malloc(3 * sizeof(double));
       
			chPerm = (int *)malloc(sizeof(int) * m->_size);
			chMinType = CS;						
			opOrder = 2;
			csmOperation(m, outAtoms, perm, &csm, dir, &dMin, CS);			

			if (csm > MINDOUBLE) {							
				for (i = 2; i <= SN_MAX; i+=2) {
					opOrder = i;
					csmOperation(m, chOutAtoms, chPerm, &chCsm, chDir, &chdMin, SN);
					if (chCsm < csm) {
						int j;
						chMinType = SN;
						csm = chCsm;
						dMin = chdMin;
						memcpy(dir, chDir, sizeof(double)* 3);
						memcpy(perm, chPerm, sizeof(int) * m->_size);	
						for (j = 0; j < m->_size; j++) { 	
							memcpy(outAtoms[j], chOutAtoms[j], sizeof(double)*3);
						}
					}
					if (csm < MINDOUBLE) break;
				}
			}
			// housekeeping
			for (i=0;i<m->_size;i++){
				free(chOutAtoms[i]);
			}
			free(chOutAtoms);	
			free(chPerm);
		}		
	}

	// print results to outFile
	if (writePDB)
		printOutputPDB(m, outAtoms, csm, dir, dMin, outFile);
	else
		printOutput(m, outAtoms, csm, dir, dMin, outFile);

	if (type == CH) {
		if (chMinType == CS) { 		
			fprintf(outFile, "\n MINIMUM CHIRALITY WAS FOUND IN CS\n\n");
		} else { 
			fprintf(outFile, "\n MINIMUM CHIRALITY WAS FOUND IN S%d\n\n", opOrder);
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

/** 
 * Compute the part of the A Matrix relevant to the current permutation 
 */
void computeMatrix(double **parray, int *perm, int size, double (*coef)[3][3], double multiplier) {
	int atom;
	for(atom = 0; atom < size; ++atom) {
		(*coef)[0][0] +=  2.0 * parray[atom][0] * parray[perm[atom]][0] * multiplier;		
		(*coef)[1][1] +=  2.0 * parray[atom][1] * parray[perm[atom]][1] * multiplier;
		(*coef)[2][2] +=  2.0 * parray[atom][2] * parray[perm[atom]][2] * multiplier;
		(*coef)[1][0] += (parray[atom][0] * parray[perm[atom]][1]+
			parray[atom][1] * parray[perm[atom]][0]) * multiplier;
		(*coef)[2][0] += (parray[atom][0] * parray[perm[atom]][2]+
			parray[atom][2] * parray[perm[atom]][0]) * multiplier;
		(*coef)[2][1] += (parray[atom][1] * parray[perm[atom]][2]+
			parray[atom][2] * parray[perm[atom]][1]) * multiplier;
	}
	(*coef)[0][1] = (*coef)[1][0];
	(*coef)[0][2] = (*coef)[2][0];
	(*coef)[1][2] = (*coef)[2][1];
}


/** 
* Compute the part of the B vector relevant to the current permutation 
*/
void computeVector(double **parray, int *perm, int size, double (*vector)[3], double multiplier) {
	int atom;
	for(atom = 0; atom < size; ++atom) {
		(*vector)[0] += multiplier * (parray[atom][1] * parray[perm[atom]][2] - parray[atom][2] * parray[perm[atom]][1]);
		(*vector)[1] += multiplier * (parray[atom][2] * parray[perm[atom]][0] - parray[atom][0] * parray[perm[atom]][2]);
		(*vector)[2] += multiplier * (parray[atom][0] * parray[perm[atom]][1] - parray[atom][1] * parray[perm[atom]][0]);
	}
}

/**
 * Calculate the best axis, and compute the CSM for it, given a pairing of the indices (perm)
 */ 
double calcRefPlane(Molecule* m, int* perm, double *dir, OperationType type) {
	double **copyMat = dmatrix(1,3,1,3);
	double *copyVec = dvector(1,3);
	double *diag = dvector(1,3);
	double *secdiag = dvector(1,3);
	double *temp = dvector(1,3);

	double matrix[3][3] = { {0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
	double vector[3] = {0.0,0.0,0.0};
	double scalar[3];
	double maxval, scl, angle;
	double rtr[6], rti[6];
	double coeffs[7];
	int i,j;
	int *curPerm = (int *)malloc(sizeof(int) * m->_size);
	double csm, dists;

	int isImproper = (type != CN) ? TRUE : FALSE;
	int isZeroAngle = (type == CS) ? TRUE : FALSE;

	// initialize identity permutation
	for (i = 0; i < m->_size; i++) {
		curPerm[i] = i;
	}

	// compute matrices according to current perm and its powers (the identity does not contribute anyway)
	for (i = 1; i < opOrder; i++) {			
		angle = isZeroAngle ? 0.0 : (2 * PI * i / opOrder);
		// The i'th power of the permutation
		for (j = 0; j < m->_size; j++) {
			curPerm[j] = perm[curPerm[j]];
		}
		if (isImproper && ((i % 2) == 1)) {
			computeMatrix(m->_pos, curPerm, m->_size, &matrix, -1-cos(angle));	
		} else {
			computeMatrix(m->_pos, curPerm, m->_size, &matrix, 1-cos(angle));	
		}
		computeVector(m->_pos, curPerm, m->_size, &vector, sin(angle));
	}
	
	// perhaps actual copying is needed
	for (i = 0; i < 3; i++) {
		copyVec[i + 1] = vector[i];
		for (j = 0; j < 3; j++) {
			copyMat[i + 1][j + 1] = matrix[i][j];
		}
	}		

	// compute the matrix's eigenvalues and eigenvectors.
	tred2(copyMat, 3, diag, secdiag);
	tqli(diag, secdiag, 3, copyMat);

	// compute square of scalar multiplications of eigen vectors with b
	for (i = 0; i < 3; i++) {
		scalar[i] = 0.0;
		for (j = 0; j < 3; j++) {
			scalar[i] += copyVec[j + 1] * copyMat[j+1][i+1];
		}
		temp[i + 1] = scalar[i] * scalar[i];
	}

	// build the polynomial
	coeffs[0] = 1.0;	// x^6
	coeffs[1] = -2 * (diag[1] + diag[2] + diag[3]);	// x^5
	coeffs[2] = diag[1] * diag[1] + diag[2] * diag[2] + diag[3] * diag[3] - 
		temp[1] - temp[2] - temp[3] + 
		4 * (diag[1] * diag[2] + diag[1] * diag[3] + diag[2] * diag[3]); // x^4
	coeffs[3] = -8 * diag[1] * diag[2] * diag[3] + 
		2 * (temp[1] * diag[2] + 
		temp[1] * diag[3] + 
		temp[2] * diag[1] + 
		temp[2] * diag[3] + 
		temp[3] * diag[1] + 
		temp[3] * diag[2] - 
		diag[1] * diag[3] * diag[3] -
		diag[1] * diag[1] * diag[2] -
		diag[1] * diag[1] * diag[3] - 
		diag[1] * diag[2] * diag[2] -
		diag[2] * diag[2] * diag[3] - 
		diag[2] * diag[3] * diag[3]); // x^3
	coeffs[4] = 4 * 
		((diag[1] * diag[2] * diag[3] * (diag[1] + diag[2] + diag[3]) -
		(temp[3] * diag[1] * diag[2] + 
		temp[2] * diag[1] * diag[3] +
		temp[1] * diag[3] * diag[2]))) - 
		temp[1] * (diag[2] * diag[2] + diag[3] * diag[3]) - 
		temp[2] * (diag[1] * diag[1] + diag[3] * diag[3]) - 
		temp[3] * (diag[1] * diag[1] + diag[2] * diag[2]) + 
		diag[1] * diag[1] * diag[2] * diag[2] +   
		diag[2] * diag[2] * diag[3] * diag[3] + 
		diag[1] * diag[1] * diag[3] * diag[3]; // x^2
	coeffs[5] = 2 * 
		(temp[1] * diag[2] * diag[3] * (diag[2] + diag[3]) + 
		temp[2] * diag[1] * diag[3] * (diag[1] + diag[3]) + 
		temp[3] * diag[1] * diag[2] * (diag[1] + diag[2])) 
		- 2 * 
		(diag[1] * diag[2] * diag[2] * diag[3] * diag[3] + 
		diag[1] * diag[1] * diag[2] * diag[3] * diag[3] + 
		diag[1] * diag[1] * diag[2] * diag[2] * diag[3]); // x
	coeffs[6] = - temp[1] * diag[2] * diag[2] * diag[3] * diag[3] - 
		temp[2] * diag[1] * diag[1] * diag[3] * diag[3] - 
		temp[3] * diag[1] * diag[1] * diag[2] * diag[2] + 
		diag[1] * diag[1] * diag[2] * diag[2] * diag[3] * diag[3]; // 1

	// solve polynomial and find maximum eigenvalue and eigen vector
	rpoly(coeffs, 6, rtr, rti);

	maxval = -MAXDOUBLE;
	for (i = 0; i < 6; i++) {				
		if (maxval < rtr[i] && (fabs(rti[i]) < 1e-6)) maxval = rtr[i];
	}			

	scl = 0.0;
	for (i = 1; i <= 3; i++) {	
		dir[i - 1]= 0.0;
		for (j = 1; j <=3; j++) {
			if ((fabs(diag[j] - maxval) < MINDOUBLE)) {
				dir[i - 1]= copyMat[i][j];
				break;
			} else {
				
				dir[i - 1] += scalar[j - 1] / (diag[j] - maxval) * copyMat[i][j];
			}
		}		
				
		scl += dir[i - 1] * copyVec[i];
	}	

	// initialize identity permutation
	for (i = 0; i < m->_size; i++) {
		curPerm[i] = i;
	}

	// compute CSM	
	csm = 0.0;		
	for (i = 0; i < opOrder; i++) {		
		// This can be more efficient - save results of matrix computation
		if (i != 0) {
			angle = isZeroAngle ? 0.0 : (2 * PI * i / opOrder);
			dists = 0.0;
			for (j = 0; j < m->_size; j++) {
				// i'th power of permutation
				curPerm[j] = perm[curPerm[j]];
			}
			for (j = 0; j < m->_size; j++) {
				dists += (m->_pos[j][0] * m->_pos[curPerm[j]][0] + 
					m->_pos[j][1] * m->_pos[curPerm[j]][1] + 
					m->_pos[j][2] * m->_pos[curPerm[j]][2]); 
			}
			csm += cos(angle) * dists;

		} else {
			csm += 1.0;	
		}
	}	

	csm += (maxval -  scl) / 2;
	csm = fabs(100 * (1.0 - csm / opOrder));
	free(curPerm);
	free_dmatrix(copyMat, 1, 3, 1, 3);
	free_dvector(copyVec, 1, 3);
	free_dvector(diag, 1, 3);
	free_dvector(secdiag, 1, 3);
	free_dvector(temp,1,3);
	return csm;
}

double createSymmetricStructure(Molecule* m, double **outAtoms, int *perm, double *dir, OperationType type, double dMin) {
	int isImproper = (type != CN) ? TRUE : FALSE;
	int isZeroAngle = (type == CS) ? TRUE : FALSE;
	int i, j, k, l;
	int *curPerm = (int *)malloc(sizeof(int) * m->_size);
	double rotaionMatrix[3][3];
	double tmpMatrix[3][3] = {{0.0, -dir[2], dir[1]}, {dir[2], 0.0, -dir[0]}, {-dir[1], dir[0], 0.0}};
	double angle;
	double res = 0.0;

	for (i = 0; i < m->_size; i++) {
		curPerm[i] = i;
		for (j = 0; j < 3; j++) {
			// initialize with identity operation
			outAtoms[i][j] = m->_pos[i][j];
		}
	}

	for (i = 1; i < opOrder; i++) {
		angle = isZeroAngle ? 0.0 : (2 * PI * i / opOrder);
		int factor = ((isImproper && (i % 2) == 1) ? (-1) : 1);
		for (j = 0; j < m->_size; j++) {
			curPerm[j] = perm[curPerm[j]];
		}
		for (j = 0; j < 3; j++) {
			for (k = 0; k < 3; k++) {
				rotaionMatrix[j][k] = 
					((j == k) ? cos(angle) : 0) + 
					(factor - cos(angle)) * dir[j] * dir[k] + 
					sin(angle) * tmpMatrix[j][k];
			}
		}
		for (j = 0; j < m->_size; j++) {
			for (k = 0; k < 3; k++) {
				for (l = 0; l < 3; l++) {
					outAtoms[curPerm[j]][k] += rotaionMatrix[k][l] * m->_pos[j][l];
				}
			}
		}
	}
	for (j = 0; j < m->_size; j++) {
		for (k = 0; k < 3; k++) {
			outAtoms[j][k] /= opOrder;
			outAtoms[j][k] *= dMin;
			res += SQR(outAtoms[j][k]);
		}
	}
	free(curPerm);

	return sqrt(res);
}


/*
* Calculates minimal csm, dMin and directional cosines by applying the chiralityFunction
* breaks down Molecule into groups of similar atoms and calculates the above only for
* permutations that keep the similar atoms within the group ( groupPermuter class )
* once it finds the optimal permutation , calls the chiralityFunction on the optimal permutation
*/
void csmOperation(Molecule* m, double** outAtoms, int *optimalPerm, double* csm, double* dir, double* dMin, OperationType type){

	int i;
	double curCsm;
	double curDir[3];
	int *idxToPos, *posToIdx;
	int * groupSizes;
	groupPermuter* gp;

	// These are the 
	int *realPerm = (int *)malloc(sizeof(int) * m->_size);

	//normalize Molecule
	if (!normalizeMolecule(m)){
		printf("Failed to normalize atom positions: dimension of set of points = zero\n");
		exit(1);
	}

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

	// create permuter
	gp = createGroupPermuter(m->_groupNum,groupSizes,m->_size,opOrder);
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

	//normalize Molecule
	if (!normalizeMolecule(m)){
		if (writeOpenu) {
			printf("ERR* Failed to normalize atom positions: dimension of set of points = zero *ERR\n");
		} else {
			printf("Failed to normalize atom positions: dimension of set of points = zero\n");
		}
		exit(1);
	}

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
	printf("%s: %.4lf\n",opName,fabs(csm));
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

	fprintf(out, "\n DIRECTIONAL COSINUSES:\n\n");
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
	printf( "DIRECTIONAL COSINUSES: %lf %lf %lf\n", dir[0], dir[1], dir[2]);

}



