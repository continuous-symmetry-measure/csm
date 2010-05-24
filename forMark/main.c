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
void csmOperation(Molecule* m, double** outAtoms, int *optimalPerm, double* csm, double* dir, double* dMin, OperationType type);
void initIndexArrays(Molecule* m, int* posToIdx, int* idxToPos);
double createSymmetricStructure(Molecule* m, double **outAtom, int *perm, double *dir, OperationType type, double dMin);

// global options
int ignoreHy = FALSE;
int removeHy = FALSE;
int ignoreSym = FALSE;
int readPDB = FALSE;
int writeOpenu = FALSE;
OperationType type;
int opOrder;

// file pointers
FILE* inFile = NULL;

char opName[100];

void usage(char *op) {
	printf("Usage: %s n input_file [-options]\n", op);	
	printf("Available options are:\n");
	printf("-ignoreHy - ignore Hydrogen atoms in computations\n");
	printf("-removeHy - remove Hydrogen atoms in computations,\n");
	printf("	rebuild molecule without them and compute\n");
	printf("-ignoreSym - Ignore all atomic symbols,\n");
	printf("	performing a purely geometric operation\n");
	printf("-readPDB - Use PDB as the input format, instead of XYZ\n");
	printf("-writeOpenu - Write output in open university format\n");	
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
	return numPermutations(m, opOrder);	
}

/*
* parses the command line parameters
*/
void parseInput(int argc, char *argv[]){

	// get commandline flags
	int i;

	// try to open infile for reading
	char* inFileName = argv[2];

	// check number of arguments
	if (argc < 3){
		usage(argv[0]);
		exit(1);
	}
	opOrder = atoi(argv[1]);
		
	if ((inFile = fopen(inFileName, "rt")) == NULL){
		if (writeOpenu) {
			printf("ERR* Failed to open data file %s *ERR\n", inFileName);
		} else {
			printf("Failed to open data file %s\n", inFileName);
		}
		exit(1);
	}



	for ( i=3;  i< argc ;  i++ ){
		if (strcmp(argv[i],"-ignoreHy" ) == 0 )
			ignoreHy = TRUE;
		else if (strcmp(argv[i],"-removeHy" ) == 0 )
			removeHy = TRUE;
		else if (strcmp(argv[i],"-ignoreSym" ) == 0 )
			ignoreSym = TRUE;
		else if (strcmp(argv[i],"-readPDB" ) == 0 )
			readPDB = TRUE;
		else if (strcmp(argv[i], "-help") == 0) {
			usage(argv[0]);
			exit(0);
		}
	}
	if (writeOpenu) {
		readPDB = TRUE;		
	}
}

// ************************************************************
//       main program
// ************************************************************

/*
* main funciton - check valid parameters, parse molecule and call chirality Operation
*/
int main(int argc, char *argv[]){

	int i,j;
	double csm, dMin;
	double **outAtoms;                 // output atom coordinates
	double dir_cn[3] = {0.0, 0.0, 0.0};   // directional cosines
	double dir_cs[3] = {0.0, 0.0, 0.0};   // directional cosines
	int *perm_cn = NULL;	
	int *perm_cs = NULL;
	int **perm_matrix = NULL;
	int nSize = 0;

	// try to read molecule from infile
	Molecule* m;

	// init options
	parseInput(argc,argv);

	
	if (readPDB)
    	m = createMoleculePDB(inFile,stdout,ignoreSym);
	else
    	m = createMolecule(inFile,stdout,ignoreSym);
	if (!m){
		if (writeOpenu) {
			printf("ERR* Failed to read molecule from data file *ERR\n");
		} else {
			printf("Failed to read molecule from data file \n");
		}
		exit(1);
	}

	// strip unwanted atoms if needbe
    	if ((ignoreHy || removeHy)){
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


	// allocate memory for outAtoms
	outAtoms = (double **)malloc(m->_size * sizeof(double*));
	for (i=0;i<m->_size;i++)
		outAtoms[i] = (double *)malloc(3 * sizeof(double));

	// perform operation
	perm_cn = (int *)malloc(sizeof(int) * m->_size);
	perm_cs = (int *)malloc(sizeof(int) * m->_size);
	nSize = opOrder;	
	
	// perform operation	
	csmOperation(m, outAtoms, perm_cn, &csm, dir_cn, &dMin, CN);			
	opOrder = 2;
	csmOperation(m, outAtoms, perm_cs, &csm, dir_cs, &dMin, CS);


	// Print result
	printf("The permutation for CN is:");
	for (i = 0; i < m->_size; i++) {
		printf("%d ", perm_cn[i] + 1);
	}
	printf("\n");
	printf("And the axis is: (%4.2f, %4.2f, %4.2f)\n", dir_cn[0],dir_cn[1],dir_cn[2]);		

	printf("The permutation for CS is:");
	for (i = 0; i < m->_size; i++) {
		printf("%d ", perm_cs[i] + 1);
	}
	printf("\n");
	printf("And the axis is: (%4.2f, %4.2f, %4.2f)\n", dir_cs[0],dir_cs[1],dir_cs[2]);		

	printf("The cosine of the angle is: %4.2f\n", dir_cs[0] * dir_cn[0] + dir_cs[1] * dir_cn[1] +dir_cs[2] * dir_cn[2]);

	// Create the permutation matrix
	perm_matrix = (int **)malloc(sizeof(int *) * nSize * 2);	
	for (i = 0; i < nSize * 2; i++) {
		perm_matrix[i] = (int *)malloc(sizeof(int) * m->_size);
	}
	
	for (i = 0; i < m->_size; i++) {
		perm_matrix[0][i] = i;
		perm_matrix[nSize][i] = perm_cs[i];
	}

	for (i = 1; i < nSize; i++) {
		for (j = 0; j < m->_size; j++) { 	
			perm_matrix[i][j] = perm_matrix[i-1][perm_cn[j]];
			perm_matrix[nSize + i][j] = perm_cn[perm_matrix[nSize + i - 1][j]];
		}
	}
	
	printf("The matrix: \n");	
	for (i = 0; i < nSize * 2; i++) {		
		if (i == 0) {
			printf("E\t\t\t");
		} else if (i < nSize) { 
			printf("(C%d)^%d\t\t\t", nSize, i);
		} else if (i == nSize) {
			printf("Cs\t\t\t"); 
		} else {
			printf("Cs*(C%d)^%d\t\t", nSize, i - nSize);
		}
		for (j = 0; j < m->_size; j++) { 	
			printf("%d ",perm_matrix[i][j] + 1);
		}
		printf("\n");
	}
	
		
	// housekeeping
	for (i=0;i<m->_size;i++){
		free(outAtoms[i]);
	}

	for (i = 0; i<nSize*2; i++) 
		free(perm_matrix[i]);

	free(perm_matrix);
	free(outAtoms);
	freeMolecule(m);
	free(perm_cs);
	free(perm_cn);

	fclose(inFile);	

	return 0;
}

/*
* creates position to index, and index to position translation arrays
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
		int factor = ((isImproper && (i % 2) == 1) ? (-1) : 1);
		angle = isZeroAngle ? 0.0 : (2 * PI * i / opOrder);
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



