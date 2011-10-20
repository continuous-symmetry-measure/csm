/*
 * Author: Amir Zait
 *
 * This is the main program. 
 * It includes the algorithm itself, as well as dealing with the program's IO
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
#include <vector>

#include "f2c.h"                     /*   Mark  */

extern "C" {
int icosahedron_(integer *, integer *, double *, double *, double *);
}

#define CSMFORMAT "CSM"
#define MAXDOUBLE  100000000.0
#define MINDOUBLE  1e-8
#define GROUPSIZE_LIMIT 15
#define GROUPSIZE_FACTOR 1.32e11
#define APPROX_RUN_PER_SEC 8e4
#define TRUE 1
#define FALSE 0
#define ZERO_IM_PART_MAX (1e-3)
#define MIN_GROUPS_FOR_OUTLIERS 10

/* borrowed from libc/misc/drand48.c in Linux libc-5.4.46 this quick
* hack by Martin Hamilton <martinh@gnu.org> to make Squid build on
* Win32 with GNU-Win32 - sorry, folks! */

#ifdef _WIN32

#define N	16
#define MASK	((unsigned)(1 << (N - 1)) + (1 << (N - 1)) - 1)
#define LOW(x)	((unsigned)(x) & MASK)
#define HIGH(x)	LOW((x) >> N)
#define MUL(x, y, z)	{ long l = (long)(x) * (long)(y); \
	(z)[0] = LOW(l); (z)[1] = HIGH(l); }
#define CARRY(x, y)	((long)(x) + (long)(y) > MASK)
#define ADDEQU(x, y, z)	(z = CARRY(x, (y)), x = LOW(x + (y)))
#define X0	0x330E
#define X1	0xABCD
#define X2	0x1234
#define A0	0xE66D
#define A1	0xDEEC
#define A2	0x5
#define C	0xB

static void next(void);
static unsigned x[3] =
{X0, X1, X2}, a[3] =
{A0, A1, A2}, c = C;

double drand48(void);

double
drand48(void)
{
	static double two16m = 1.0 / (1L << N);
	next();
	return (two16m * (two16m * (two16m * x[0] + x[1]) + x[2]));
}

static void
next(void)
{
	unsigned p[2], q[2], r[2], carry0, carry1;

	MUL(a[0], x[0], p);
	ADDEQU(p[0], c, carry0);
	ADDEQU(p[1], carry0, carry1);
	MUL(a[0], x[1], q);
	ADDEQU(p[1], q[0], carry0);
	MUL(a[1], x[0], r);
	x[2] = LOW(carry0 + carry1 + CARRY(p[1], r[0]) + q[1] + r[1] +
		a[0] * x[2] + a[1] * x[1] + a[2] * x[0]);
	x[1] = LOW(p[1] + r[0]);
	x[0] = LOW(p[0]);
}

#endif /* HAVE_DRAND48 */


typedef enum {
	CN,
	SN, 
	CS,
	CI, 
	CH
} OperationType;

struct distRecord {
	double distance;
	int row; 
	int col;
};

int distComp(const void *pd1, const void* pd2) { 
	const struct distRecord* d1 = (const struct distRecord*)pd1;
	const struct distRecord* d2 = (const struct distRecord*)pd2;
	if (d1->distance > d2->distance) { 
		return 1;
	} else if (d1->distance < d2->distance) { 
		return -1;
	}
	return 0;
}

extern "C" {

// from tqli.c - nrbook
double tqli(double d[], double e[], int n, double **z);
// from tred2.c - nrbook
void tred2(double **a, int n, double d[], double e[]);
// from rpoly.c Jenkins-Traub Polynomial Solver
int rpoly(double *op, int degree, double *zeror, double *zeroi);

}

// function declarations
void printOutput(Molecule* m, double** outAtoms, double csm, double *dir, double dMin, FILE *out, double* localCSM);
void printOutputPDB(Molecule* m, double** outAtoms, double csm, double *dir, double dMin, FILE *out);
void printOutputFormat(Molecule* m, OBMol& mol, double** outAtoms, double csm, double *dir, double dMin, FILE *out, char *fname, double* localCSM);
void csmOperation(Molecule* m, double** outAtoms, int *optimalPerm, double* csm, double* dir, double* dMin, OperationType type);
void runSinglePerm(Molecule* m, double** outAtoms, int* perm, double* csm, double* dir, double* dMin, OperationType type);
void findBestPerm(Molecule* m, double** outAtoms, int* optimalPerm, double* csm, double* dir, double* dMin, OperationType type);
void findSymmetryDirection(Molecule *m, double  ***dirs, int *n_dirs, OperationType type);
void estimatePerm(Molecule* m, int *perm, double *dir, OperationType type);
void readPerm(FILE* permfile, int* perm, int size);
void lineFit(double **points, int nPoints, double **dirs, int* outlies);
void planeFit(double **points, int nPoints, double **dirs, int *outliers);
void initIndexArrays(Molecule* m, int* posToIdx, int* idxToPos);
double createSymmetricStructure(Molecule* m, double **outAtom, int *perm, double *dir, OperationType type, double dMin);
double computeLocalCSM(Molecule* m, double *localCSM, int *perm, double *dir, OperationType type);
double findSecondaryPerm(Molecule* m, double** outAtoms, int *optimalPerm, double* dir, 
			double* dMin, OperationType type, double* dir_cn, int skipIdentity);
void findBestSecondaryPerm(Molecule* m, double** outAtoms, int* perm, double* csm, 
						   double* dir, double* dMin, OperationType type, double *dir_cn, int skipIdentity);

/*------    Mark -------------------*/
int normalizeMolecule2(Molecule *m);


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
int findPerm = FALSE; 
int useMass    = FALSE;
int limitRun = TRUE;
char *format = NULL;
int babelBond = FALSE;
int timeOnly = FALSE;
int sn_max = 4;
//int anneal = FALSE;
int detectOutliers = FALSE;
double A = 2;
int babelTest = FALSE;
int printNorm = FALSE;
int printLocal = FALSE;

// Dn / CnV stuff
#ifdef CNV
OperationType secondOpType = CS;
#else // DN
OperationType secondOpType = CN;
#endif 

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
	printf("-findperm	   - Attempt to search for a permutation\n");
	printf("-detectOutliers - Use statistical methods to try and improve -findperm's results\n");
//	printf("-anneal		   - Try to anneal the result\n");
	printf("-babelbond	   - Let openbabel compute bonding\n");
	printf("-useMass	   - Use the atomic masses to define center of mass\n");
	printf("-timeOnly	   - Only print the time and exit\n");
	printf("-babelTest	   - Test if the molecule is legal or not\n");
	printf("-approx			- Equivalent to -detectOutliers -findperm together");
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
	if (operationOrder > 2 && t == SN) {	
		// In this case - we enumerate over groups of 1, 2, N
		for (i = 1; i <= m->_groupNum; i++) {			
			int groupSize = getGroupSize(m, i);
			double temp = 0;
			double fact = factorial(groupSize);
			// Enumerate over the number of groups of size two
			for (k = 0; k <= (groupSize / 2); k++) {							
				// Enumerate over the number of groups of size operationOrder
				for (j = 0; j <= ((groupSize - 2 * k) / operationOrder); j++) {				
					temp += fact/(factorial(groupSize - j * operationOrder - k * 2) *
							pow((double)operationOrder, j) *factorial(j) * 
							pow((double)2, k) * factorial(k));
				}
			}
			total *= (temp);
		}		
	} else { 
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
	}
	return total;
}

double totalNumPermutations(Molecule *m) {
	if (type != CH) { 
		return numPermutations(m, opOrder, type);	
	} else {
		// CS
		int i;
		double numPerms = numPermutations(m, 2, SN);
		for (i = 2; i <= sn_max; i+=2) {
			numPerms += numPermutations(m, i, SN);
		}	
		return numPerms;		
	}
}

char *getExtension(char *fname) {
	return strrchr(fname,'.') + 1;
}

// Math utils

/*
 * Normalizes the position of atoms of the molecule
 * returns one [TRUE] if successful, zero[FALSE] otherwise
 */
void normalize(double **coords, Molecule *m){

	double tmp,x_avg, y_avg, z_avg,norm;
	int i;

	x_avg = y_avg = z_avg = 0.0;

	double mass_sum = 0;
	for(i=0; i< m->_size; i++){
		x_avg += coords[i][0] * m->_mass[i];
		y_avg += coords[i][1] * m->_mass[i];
		z_avg += coords[i][2] * m->_mass[i];
		mass_sum += m->_mass[i];
	}
	x_avg /= (double)(mass_sum);
	y_avg /= (double)(mass_sum);
	z_avg /= (double)(mass_sum);

	norm = 0.0;
	for(i=0; i< m->_size; i++){
		tmp = SQR(coords[i][0]-x_avg) +
		      SQR(coords[i][1]-y_avg) +
		      SQR(coords[i][2]-z_avg);
		norm += tmp;
	}
	// normalize to 1 and not molecule size
	//norm = sqrt(norm / (double)m->_size);
	norm = sqrt(norm);


	for(i=0; i< m->_size; i++){
		coords[i][0] = ((coords[i][0] - x_avg) / norm);
		coords[i][1] = ((coords[i][1] - y_avg) / norm);
		coords[i][2] = ((coords[i][2] - z_avg) / norm);
	}
}

double Magnitude( double *Point1, double *Point2 )
{
	double Vector[3];

	Vector[0] = Point2[0] - Point1[0];
	Vector[1] = Point2[1] - Point1[1];
	Vector[2] = Point2[2] - Point1[2];

	return sqrt( Vector[0] * Vector[0] + Vector[1]* Vector[1] + Vector[2] * Vector[2] );
}

double computeDistanceFromLine( double *Point, double *LineStart, double *LineEnd)
{
	double LineMag;
	double U;
	double Intersection[3];

	LineMag = Magnitude( LineEnd, LineStart );

	U = ( ( ( Point[0] - LineStart[0]) * ( LineEnd[0] - LineStart[0] ) ) +
		( ( Point[1] - LineStart[1] ) * ( LineEnd[1] - LineStart[1] ) ) +
		( ( Point[2] - LineStart[2] ) * ( LineEnd[2] - LineStart[2] ) ) ) /
		( LineMag * LineMag );

	Intersection[0] = LineStart[0] + U * ( LineEnd[0] - LineStart[0] );
	Intersection[1]= LineStart[1] + U * ( LineEnd[1] - LineStart[1] );
	Intersection[2] = LineStart[2] + U * ( LineEnd[2] - LineStart[2] );

	return Magnitude( Point, Intersection );	
}


/*
*  This Quickselect routine is based on the algorithm described in
*  "Numerical recipes in C", Second Edition,
*  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
*  This code by Nicolas Devillard - 1998. Public domain.
*/


#define ELEM_SWAP(a,b) { register double t=(a);(a)=(b);(b)=t; }

double findMedian(double arr[], int n) 
{
	int low, high ;
	int median;
	int middle, ll, hh;

	low = 0 ; high = n-1 ; median = (low + high) / 2;
	for (;;) {
		if (high <= low) /* One element only */
			return arr[median] ;

		if (high == low + 1) {  /* Two elements only */
			if (arr[low] > arr[high])
				ELEM_SWAP(arr[low], arr[high]) ;
			return arr[median] ;
		}

		/* Find median of low, middle and high items; swap into position low */
		middle = (low + high) / 2;
		if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]) ;
		if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]) ;
		if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]) ;

		/* Swap low item (now in position middle) into position (low+1) */
		ELEM_SWAP(arr[middle], arr[low+1]) ;

		/* Nibble from each end towards middle, swapping items when stuck */
		ll = low + 1;
		hh = high;
		for (;;) {
			do ll++; while (arr[low] > arr[ll]) ;
			do hh--; while (arr[hh]  > arr[low]) ;

			if (hh < ll)
				break;

			ELEM_SWAP(arr[ll], arr[hh]) ;
		}

		/* Swap middle item (in position low) back into correct position */
		ELEM_SWAP(arr[low], arr[hh]) ;

		/* Re-set active partition */
		if (hh <= median)
			low = ll;
		if (hh >= median)
			high = hh - 1;
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

	opOrder = atoi(argv[1]);


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
		} else if (nextIsMaxSn) { 
			sn_max = atoi(argv[i]);
			nextIsMaxSn = FALSE;
		} else if (strcmp(argv[i],"-sn_max" ) == 0) {
			if (type != CH) { 
				printf("This option only applies to chirality\n");
				exit(1);
			}
			nextIsMaxSn = TRUE;	
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
		} else if  (strcmp(argv[i], "-printNorm") == 0) {
			printNorm = TRUE;
		} else if (strcmp(argv[i], "-help") == 0) {
			usage(argv[0]);
			exit(0);
		} else if (strcmp(argv[i], "-findperm") == 0) { 
			findPerm = TRUE;
		} else if (strcmp(argv[i], "-detectOutliers") == 0) {
			detectOutliers = TRUE;
		} else if (strcmp(argv[i], "-approx") == 0) {
			detectOutliers = TRUE;
			findPerm = TRUE;
		} else if (strcmp(argv[i], "-babelTest") == 0) { 
			babelTest = TRUE;
		} else if (strcmp(argv[i], "-printlocal") == 0) { 
			printLocal = TRUE;
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

	
	int **perm_matrix = NULL;
	int i,j,nSize;
	double csm,csm_cn,csm_sec, dMin;
	double **outAtoms;                 // output atom coordinates
	double dir_cn[3] = {0.0, 0.0, 0.0};   // directional cosines
	double dir_sec[3] = {0.0, 0.0, 0.0};   // directional cosines
	int *perm_cn = NULL;	
	int *perm_sec = NULL;
	double *localCSM = NULL;		 
	int skipIdentity = 0;

	// init options
	parseInput(argc,argv);

	if (findPerm && useperm) { 
		printf("-findperm and -useperm can't be used together...");
		exit(1);
	} 

	// try to read molecule from infile
	Molecule* m;
	OBMol mol; 
	
	if (useFormat) {
		// If a specific format is used, read molecule using that format
		if (strcasecmp(format, CSMFORMAT) == 0) {
			m = createMolecule(inFile,stdout,ignoreSym && !useperm);
			if (m==NULL) exit(1);
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
			if (m==NULL) exit(1);
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

	if (babelTest) // Mol is ok - return 0
		return 0;

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

	// allocate memory for outAtoms
	outAtoms = (double **)malloc(m->_size * sizeof(double*));
	for (i=0;i<m->_size;i++)
		outAtoms[i] = (double *)malloc(3 * sizeof(double));
       

	//normalize Molecule
	if (!normalizeMolecule(m,FALSE)){
		if (writeOpenu) {
			printf("ERR* Failed to normalize atom positions: dimension of set of points = zero *ERR\n");
		} else {
			printf("Failed to normalize atom positions: dimension of set of points = zero\n");
		}
		exit(1);
	}

	// perform operation
	perm_cn = (int *)malloc(sizeof(int) * m->_size);
	perm_sec = (int *)malloc(sizeof(int) * m->_size);
	nSize = opOrder;	
	skipIdentity = (secondOpType == CS) && (opOrder > 2);
	
	// perform operation	
	normalize(outAtoms, m);	
	if (findPerm) {
		findBestPerm(m, outAtoms, perm_cn, &csm, dir_cn, &dMin, CN);			
	} else {
		csmOperation(m, outAtoms, perm_cn, &csm, dir_cn, &dMin, CN);			
	}
	csm_cn = csm;
	opOrder = 2;
	if (findPerm) {
	    findBestSecondaryPerm(m, outAtoms, perm_sec, &csm_sec, dir_sec, &dMin, secondOpType, dir_cn,skipIdentity);
	} else {
		csm_sec = findSecondaryPerm(m, outAtoms, perm_sec, dir_sec, &dMin, secondOpType, dir_cn, skipIdentity);
	}

	// Print result
	
	fprintf(outFile,"The permutation for CN is:");
	for (i = 0; i < m->_size; i++) {
		fprintf(outFile,"%d ", perm_cn[i] + 1);
	}
	fprintf(outFile,"\n");
	fprintf(outFile,"And the axis is: (%4.2f, %4.2f, %4.2f)\n", dir_cn[0],dir_cn[1],dir_cn[2]);		
	
	printf("The CN axis is: (%4.2f, %4.2f, %4.2f)\n", dir_cn[0],dir_cn[1],dir_cn[2]);		

#ifdef CNV
	fprintf(outFile,"The permutation for CS is:");
#else // DN
	fprintf(outFile,"The permutation for C2 is:");
#endif
	for (i = 0; i < m->_size; i++) {
		fprintf(outFile,"%d ", perm_sec[i] + 1);
	}
	fprintf(outFile,"\n");
	fprintf(outFile,"And the axis is: (%4.2f, %4.2f, %4.2f)\n", dir_sec[0],dir_sec[1],dir_sec[2]);		

#ifdef CNV
	printf("The CS axis is: (%4.2f, %4.2f, %4.2f)\n", dir_sec[0],dir_sec[1],dir_sec[2]);		
#else // DN
	printf("The C2 axis is: (%4.2f, %4.2f, %4.2f)\n", dir_sec[0],dir_sec[1],dir_sec[2]);		
#endif

	printf("The cosine of the angle is: %4.2f\n", dir_sec[0] * dir_cn[0] + dir_sec[1] * dir_cn[1] +dir_sec[2] * dir_cn[2]);
	fprintf(outFile,"The cosine of the angle is: %4.2f\n", dir_sec[0] * dir_cn[0] + dir_sec[1] * dir_cn[1] +dir_sec[2] * dir_cn[2]);

	printf("CSM for cn: %4.2f\n", csm_cn);
	fprintf(outFile,"CSM for cn: %4.2f\n", csm_cn);

#ifdef CNV
	printf("CSM for cs: %4.2f\n", csm_sec);
	fprintf(outFile,"CSM for cs: %4.2f\n", csm_sec);
#else // DN
	printf("CSM for c2: %4.2f\n", csm_sec);
	fprintf(outFile,"CSM for c2: %4.2f\n", csm_sec);
#endif

	// Create the permutation matrix
	perm_matrix = (int **)malloc(sizeof(int *) * nSize * 2);	
	for (i = 0; i < nSize * 2; i++) {
		perm_matrix[i] = (int *)malloc(sizeof(int) * m->_size);
	}
	
	for (i = 0; i < m->_size; i++) {
		perm_matrix[0][i] = i;
		perm_matrix[nSize][i] = perm_sec[i];
	}

	for (i = 1; i < nSize; i++) {
		for (j = 0; j < m->_size; j++) { 	
			perm_matrix[i][j] = perm_matrix[i-1][perm_cn[j]];
			perm_matrix[nSize + i][j] = perm_cn[perm_matrix[nSize + i - 1][j]];
		}
	}
	
	fprintf(outFile,"The permutations matrix: \n");	
	for (i = 0; i < nSize * 2; i++) {			
		if (i == 0) {
			fprintf(outFile,"E\t\t\t");
		} else if (i < nSize) { 
			fprintf(outFile,"(C%d)^%d\t\t\t", nSize, i);
		} else if (i == nSize) {
#ifdef CNV
			fprintf(outFile,"CS\t\t\t"); 
#else // DN
			fprintf(outFile,"C2\t\t\t"); 
#endif 
		} else {
#ifdef CNV
			fprintf(outFile,"CS*(C%d)^%d\t\t", nSize, i - nSize);
#else // DN
			fprintf(outFile,"C2*(C%d)^%d\t\t", nSize, i - nSize);
#endif
		}
		
		
		
		for (j = 0; j < m->_size; j++) { 	
			fprintf(outFile,"%d ",perm_matrix[i][j] + 1);
		}
		fprintf(outFile,"\n");
	}
	
#ifdef CNV
	/**************    Mark begin	 *******************/
{


double calc_cnv(Molecule* m, double **outAtoms, double *dir, int nSize,
              int **perm_matrix);

    
    
double *x0, *y0, *z0;
double dir[3], dir_out[3];
integer n_resolution, n_axes;
double **rotAtoms;
double sym, s = MAXDOUBLE;

//  Normalize Molecule to rms size
//  After normalization rms size of molecular is equal to 1
if (!normalizeMolecule2(m)){
   printf("Failed to normalize atom positions: dimension of set of points = zero\n");
   exit(1);
}  


// allocate memory for rotAtoms
rotAtoms = (double **)malloc(m->_size * sizeof(double*));
for (i=0;i<m->_size;i++)
    rotAtoms[i] = (double *)malloc(3 * sizeof(double));

x0 = (double*)malloc(sizeof(double)*500000);
y0 = (double*)malloc(sizeof(double)*500000);
z0 = (double*)malloc(sizeof(double)*500000);
	    
//  Array of "n_axes" axis	    
n_resolution = 110;
icosahedron_(&n_resolution, &n_axes, x0, y0, z0);
 

//  Calculation of Cnv CSM

for (i = 0; i < n_axes; i++) {

    
    dir[0] = x0[i]; dir[1] = y0[i]; dir[2] = z0[i]; 
    sym = calc_cnv(m, rotAtoms, dir, nSize, perm_matrix);
          
        
    if (sym < s) {
        s = sym;
	for (j = 0; j < 3; j++)
	    dir_out[j] = dir[j];
        }
     }
	    
/*
printf("\nCnv CSM is equal to    %20.10lf\n",s);	    
printf("\nRotation axes  %15.10lf %15.10lf %15.10lf\n",dir_out[0],dir_out[1],dir_out[2]);	    
   
*/
    
  

printf("\nC%dV: %6.4lf\n",nSize, abs(100.*s));    
fprintf(outFile,"\nC%dV: %6.4lf\n",nSize,abs(100.*s));    
    
    
free(x0);
free(y0);
free(z0);

}

/*****************     Mark End     *******************/

#else // DN
	
{

double calc_dn(Molecule* m, double **outAtoms, double *dir, int nSize,
              int **perm_matrix);

    
    
double *x0, *y0, *z0;
double dir[3], dir_out[3];
integer n_resolution, n_axes;
double **rotAtoms;
double sym, s = MAXDOUBLE;

//  Normalize Molecule to rms size
//  After normalization rms size of molecular is equal to 1
if (!normalizeMolecule2(m)){
   printf("Failed to normalize atom positions: dimension of set of points = zero\n");
   exit(1);
}  

// allocate memory for rotAtoms
rotAtoms = (double **)malloc(m->_size * sizeof(double*));
for (i=0;i<m->_size;i++)
    rotAtoms[i] = (double *)malloc(3 * sizeof(double));

x0 = (double*)malloc(sizeof(double)*500000);
y0 = (double*)malloc(sizeof(double)*500000);
z0 = (double*)malloc(sizeof(double)*500000);
	    
//  Array of "n_axes" axis	    
n_resolution = 110;
icosahedron_(&n_resolution, &n_axes, x0, y0, z0);
 

//  Calculation of Dn CSM

for (i = 0; i < n_axes; i++) {    
    
    dir[0] = x0[i]; dir[1] = y0[i]; dir[2] = z0[i];
    
/*    dir[0] =-1.0; dir[1] = 0.0; dir[2] = 0.0;*/

    sym = calc_dn(m, rotAtoms, dir, nSize, perm_matrix);  
        
    if (sym < s) {
        s = sym;
	for (j = 0; j < 3; j++)
	    dir_out[j] = dir[j];
        }
     }
	    
/*
printf("\nDn CSM is equal to    %20.10lf\n",s);

*/

printf("\nRotation axes  %15.10lf %15.10lf %15.10lf\n",dir_out[0],dir_out[1],dir_out[2]);	    
fprintf(outFile, "\nRotation axes  %15.10lf %15.10lf %15.10lf\n",dir_out[0],dir_out[1],dir_out[2]);	    

printf("\nD%d: %15.10lf\n",nSize,abs(100.*s));    
fprintf(outFile, "\nD%d:  %15.10lf\n",nSize,abs(100.*s));    
    
free(x0);
free(y0);
free(z0);
} 


#endif 

	// De-normalize
	for (i = 0; i < m->_size; i++) { 
		m->_pos[i][0] *= m->_norm;
		m->_pos[i][1] *= m->_norm;
		m->_pos[i][2] *= m->_norm;
		outAtoms[i][0] *= m->_norm;
		outAtoms[i][1] *= m->_norm;
		outAtoms[i][2] *= m->_norm;
	}	

	fprintf(outFile,"\n");

	// housekeeping
	for (i=0;i<m->_size;i++){
		free(outAtoms[i]);
	}
	free(outAtoms);
	freeMolecule(m);
	free(perm_cn);
	free(perm_sec);

	if (printLocal) free(localCSM);

	fclose(inFile);
	fclose(outFile);

	for (i = 0; i<nSize*2; i++) 
		free(perm_matrix[i]);

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
* Compute the part of the B vec relevant to the current permutation 
*/
void computeVector(double **parray, int *perm, int size, double (*vec)[3], double multiplier) {
	int atom;
	for(atom = 0; atom < size; ++atom) {		
		(*vec)[0] += multiplier * (parray[atom][1] * parray[perm[atom]][2] - parray[atom][2] * parray[perm[atom]][1]);
		(*vec)[1] += multiplier * (parray[atom][2] * parray[perm[atom]][0] - parray[atom][0] * parray[perm[atom]][2]);
		(*vec)[2] += multiplier * (parray[atom][0] * parray[perm[atom]][1] - parray[atom][1] * parray[perm[atom]][0]);
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
	double vec[3] = {0.0,0.0,0.0};
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
		computeVector(m->_pos, curPerm, m->_size, &vec, sin(angle));
	}
	
	// perhaps actual copying is needed
	for (i = 0; i < 3; i++) {
		copyVec[i + 1] = vec[i];
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

	// DEBUG printf("%.6f %.6f %.6f\n", copyVec[1], copyVec[2], copyVec[3]);

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

	// solve polynomial and find maximum eigenvalue and eigen vec
	rpoly(coeffs, 6, rtr, rti);

	maxval = -MAXDOUBLE;
	for (i = 0; i < 6; i++) {								
		if (maxval < rtr[i] && (fabs(rti[i]) < ZERO_IM_PART_MAX)) maxval = rtr[i];		
	}	

	scl = 0.0;
	if ((isZeroAngle) || (opOrder == 2)) {
		// If we are in zero angle case, we should pick the direction matching maxval
		double minDist = MAXDOUBLE;
		int minarg = 0;
		for (i = 1; i <= 3; i ++) {
			if (fabs(diag[i] - maxval) < minDist) {
				minDist = fabs(diag[i] - maxval);
				minarg = i;
			}
		}
		for (i = 0; i < 3; i++) { 
			dir[i] = copyMat[i+1][minarg];
		}
	} else {
		for (i = 1; i <= 3; i++) {	
			dir[i - 1]= 0.0;
			for (j = 1; j <=3; j++) {					
				// error safety
				if (fabs(diag[j] - maxval) < 1e-6) {				
					dir[i - 1]= copyMat[i][j];
					break;
				} else {				
					dir[i - 1] += scalar[j - 1] / (diag[j] - maxval) * copyMat[i][j];			
				}
			}		
			scl += dir[i - 1] * copyVec[i];
		}						
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
			angle = isZeroAngle ? 0.0 : ((2 * PI * i) / opOrder);
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

	
	// DEBUG printf("%.6f %.6f %.6f\n", csm, maxval, scl);
	// DEBUG printf("dir: %.6f %.6f %.6f\n", dir[0], dir[1], dir[2]);
	// DEBUG printf("eig: %.6f %.6f %.6f\n",  diag[1], diag[2], diag[3]);
	csm += (maxval -  scl) / 2;
	csm = fabs(100 * (1.0 - csm / opOrder));
	free(curPerm);

	//printf("%4.2f %4.2f %4.2f - %4.2f\n", dir[0], dir[1], dir[2], csm);
	
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
		// initialize with identity operation
		curPerm[i] = i;
		for (j = 0; j < 3; j++) {
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
					outAtoms[j][k] += rotaionMatrix[k][l] * m->_pos[curPerm[j]][l];
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

double computeLocalCSM(Molecule* m, double *localCSM, int *perm, double *dir, OperationType type) {
	int isImproper = (type != CN) ? TRUE : FALSE;
	int isZeroAngle = (type == CS) ? TRUE : FALSE;
	int i, j, k, l;
	int *curPerm = (int *)malloc(sizeof(int) * m->_size);
	double rotaionMatrix[3][3];
	double tmpMatrix[3][3] = {{0.0, -dir[2], dir[1]}, {dir[2], 0.0, -dir[0]}, {-dir[1], dir[0], 0.0}};
	double angle;
	double res = 0.0;

	double tempCSM = 0.0;


	for (i = 0; i < m->_size; i++) {
		// initialize with identity operation
		curPerm[i] = i;
		localCSM[i] = 0;
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
			double rotated[3] = {0.0,0.0,0.0};
			for (k = 0; k < 3; k++) {
				for (l = 0; l < 3; l++) {
					rotated[k] += rotaionMatrix[k][l] * m->_pos[curPerm[j]][l];
				}
			}			
			for (k = 0; k < 3; k++) {
				tempCSM += SQR(rotated[k] - m->_pos[j][k]);
				localCSM[j] += SQR(rotated[k] - m->_pos[j][k]);
			}			
		}

	}

	for (j = 0; j < m->_size; j++) {
		localCSM[j] *= 100 / (2 * opOrder);
	}
	tempCSM *= 100.0 / (2 * opOrder);	

	return tempCSM;
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

	// create permuter
	if (type == SN && opOrder > 2) {
		addGroupsOfTwo = 1;
	} else {
		addGroupsOfTwo = 0;
	}
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
* Calculates minimal csm, dMin and directional cosines by applying the chiralityFunction
* breaks down Molecule into groups of similar atoms and calculates the above only for
* permutations that keep the similar atoms within the group ( groupPermuter class )
* once it finds the optimal permutation , calls the chiralityFunction on the optimal permutation
*/
double findSecondaryPerm(Molecule* m, double** outAtoms, int *optimalPerm, double* dir, 
			double* dMin, OperationType type, double* dir_cn, int skipIdentity){

	int i;
	double bestOrth = MAXDOUBLE,curOrth,bestCsm;
	double curDir[3];
	int *idxToPos, *posToIdx;
	int * groupSizes;
	groupPermuter* gp;
	int addGroupsOfTwo;
	double norm,curCsm = 0;

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
	bestCsm = curCsm = MAXDOUBLE;

	// get group sizes
	for(i=1; i<= m->_groupNum ; i++){
		groupSizes[i-1] = getGroupSize(m,i);
	}

	// create permuter
	if (type == SN && opOrder > 2) {
		addGroupsOfTwo = 1;
	} else {
		addGroupsOfTwo = 0;
	}
	gp = createGroupPermuter(m->_groupNum,groupSizes,m->_size,opOrder, addGroupsOfTwo);
	if (!gp){
		if (writeOpenu) {
			printf("ERR* Failed to create groupPermuter *ERR\n");
		} else {
			printf("Failed to create groupPermuter \n");	
		}	
		exit(1);
	};

	/* Since we are looking for orthogonal - skip identity perm */
	if (skipIdentity) nextGroupPermutation(gp);

	// calculate csm for each valid permutation & remember minimal (in optimalAntimer)
	while ( nextGroupPermutation(gp) ) {

		for (i = 0; i < m->_size; i++) {
			realPerm[i] = idxToPos[gp->_index[posToIdx[i]]];			
		}			
		curCsm = calcRefPlane(m, realPerm, curDir, type);		
		norm = curDir[0] * curDir[0] + curDir[1] * curDir[1] + curDir[2] * curDir[2];
		if (norm < 0.5) continue;

		curOrth = abs(curDir[0] * dir_cn[0] + curDir[1] * dir_cn[1] + curDir[2] * dir_cn[2]);
		// Make all almost-orthogonal values completely orthogonal	
		// make only 100 different values...
		curOrth = floor(curOrth * 100) * 1.0 / 100;

		if(curOrth < bestOrth ||(curOrth == bestOrth && curCsm < bestCsm)) {
			bestOrth = curOrth;
			bestCsm = curCsm;
			dir[0] = curDir[0]; dir[1] = curDir[1]; dir[2] = curDir[2];
			
			for (i = 0; i < m->_size; i++) {
				optimalPerm[i] = realPerm[i];
			}
		}			
	}

	// failed to find value for any permutation
	
	if (bestCsm == MAXDOUBLE){
		if (writeOpenu) {
			printf("ERR* Failed to calculate a csm value for %s *ERR\n",opName);
		} else {
			printf("Failed to calculate a csm value for %s \n",opName);
		}
		exit(1);
	}

	// which is DMIN?
	*dMin = (1.0 - (bestCsm / 100 * opOrder / (opOrder - 1)));
	createSymmetricStructure(m, outAtoms, optimalPerm, dir, type, *dMin);		

	// housekeeping
	free(groupSizes);
	free(idxToPos);
	free(posToIdx);
	freeGroupPermuter(gp);
	free(realPerm);

	return bestCsm;
}


/*
 * Calculates csm, dMin and directional cosines for a given permutation
 */
void runSinglePerm(Molecule* m, double** outAtoms, int *perm, double* csm, double* dir, double* dMin, 
				   OperationType type){
	*csm = calcRefPlane(m, perm, dir, type);

	// which is DMIN?
	*dMin = (1.0 - (*csm / 100 * opOrder / (opOrder - 1)));
	createSymmetricStructure(m, outAtoms, perm, dir, type, *dMin);
}

int isIdentityPerm(int *perm, int size) {
	int i = 0; 
	for (i = 0; i < size; i++) {
		if (perm[i] != i) return 0;
	}
	return 1;
}

/**
 * Finds an approximate permutation which can be used in the analytical computation.
 */
void findBestSecondaryPerm(Molecule* m, double** outAtoms, int* perm, double* csm, double* dir, 
						   double* dMin, OperationType type, double *dir_cn, int skipIdentity) {

	int *temp = (int*)malloc(sizeof(int) * m->_size);
	int i = 0;	
	double** dirs;
	int n_dirs;
	int *bestPerm = (int*)malloc(sizeof(int) * m->_size);	
	double bestDir[3], bestOrth = MAXDOUBLE;
	int maxIters = 50;	
	*csm = MAXDOUBLE;

	// Find an initial approximated symmetry axis/plain
	findSymmetryDirection(m, &dirs, &n_dirs, type);

	// Find the permutation best matching this direction - there are n_dirs results to the algorithm		
	for (i = 0; i < n_dirs; i++) { 
		double dist = MAXDOUBLE, old = MAXDOUBLE, best = MAXDOUBLE;	
		double curOrth, norm;
		double tempDir[3];
		int iterNum = 1;
		estimatePerm(m, bestPerm, dirs[i], type);	
		runSinglePerm(m, outAtoms, bestPerm, &dist, tempDir, dMin, type);	
		memcpy(bestDir, tempDir, sizeof(double) * 3);
		old = MAXDOUBLE; best = dist;

		// solve analytically using this permutation, repeat until converged
		// Pick the best one
		// Stop once:
		// 1. The csm is less than 1e-4
		// 2. The difference between the old and new csm is less than 1%
		// 3. The max number of iterations has been reached
		while ((fabs(dist) > 1e-4) && (fabs(old - dist)/fabs(old) > 0.01) && (iterNum < maxIters)) {
			old = dist;
			estimatePerm(m, temp, tempDir, type);		
			if (skipIdentity && isIdentityPerm(temp,m->_size)) continue;
			runSinglePerm(m, outAtoms, temp, &dist, tempDir, dMin, type);					
			
			if (dist < best) {
				best = dist;
				memcpy(bestPerm, temp, sizeof(int) * m->_size);
				memcpy(bestDir, tempDir, sizeof(double) * 3);
			}	
		
			iterNum++;			
			//printf("Old csm: %6.4f New csm %6.4f\n", old, dist);
		}

		norm = bestDir[0] * bestDir[0] + bestDir[1] * bestDir[1] + bestDir[2] * bestDir[2];
		if (norm < 0.5) continue;

		curOrth = abs(bestDir[0] * dir_cn[0] + bestDir[1] * dir_cn[1] + bestDir[2] * dir_cn[2]);
		// Make all almost-orthogonal values completely orthogonal	
		// make only 100 different values...
		curOrth = floor(curOrth * 100) * 1.0 / 100;
		if (skipIdentity && isIdentityPerm(bestPerm,m->_size)) 
		{
			printf("Identity permutation found - not allowed in cnv for n>2\n");
			continue;
		}
		if((curOrth < bestOrth) || (curOrth == bestOrth && best < *csm)) { 
			bestOrth = curOrth;
			*csm = best;
			dir[0] = bestDir[0]; dir[1] = bestDir[1]; dir[2] = bestDir[2];
			for (int i = 0; i < m->_size; i++) {
				perm[i] = bestPerm[i];
			}
		}			

		printf("Attempt #%d: csm is %4.2f, cos(angle): %4.2f after %d iterations\n", (i+1), best,curOrth, iterNum);				
	}

	for (i = 0; i < n_dirs; i++) {
		free(dirs[i]);
	}

	free(dirs);
	free(bestPerm);	

	// run once more to get the out atoms right !
	runSinglePerm(m, outAtoms, perm, csm, dir, dMin, type);					
	
	free(temp);	
	
}


/**
 * Finds an approximate permutation which can be used in the analytical computation.
 */
void findBestPerm(Molecule* m, double** outAtoms, int* perm, double* csm, double* dir, double* dMin, OperationType type) {
	int *temp = (int*)malloc(sizeof(int) * m->_size);
	int i = 0;	
	// The algorithm aims to find the best perm which can then be used for the analytic solution	
	if ((type == CI) || (type == SN && opOrder == 2)) { 		
		// For inversion - simply find for each orbit the best matching - the closest after the operation.	
		// Do nothing - no need to find axis. yay.
		dir[0] = 1.0; dir[1] = 0.0; dir[2] = 0.0;
		estimatePerm(m, perm, dir, type);
		runSinglePerm(m, outAtoms, perm, csm, dir, dMin, type);	
	} else { 		
		double** dirs;
		int n_dirs;
		int *bestPerm = (int*)malloc(sizeof(int) * m->_size);	
		double bestDir[3];
		int maxIters = 50;	
		*csm = MAXDOUBLE;

		// Find an initial approximated symmetry axis/plain
		findSymmetryDirection(m, &dirs, &n_dirs, type);

		// Find the permutation best matching this direction - there are n_dirs results to the algorithm		
		for (i = 0; i < n_dirs; i++) { 
			double dist = MAXDOUBLE, old = MAXDOUBLE, best = MAXDOUBLE;	
			double tempDir[3];
			int iterNum = 1;
			estimatePerm(m, bestPerm, dirs[i], type);	
			runSinglePerm(m, outAtoms, bestPerm, &dist, tempDir, dMin, type);	
			memcpy(bestDir, tempDir, sizeof(double) * 3);
			old = MAXDOUBLE; best = dist;
			
			// solve analytically using this permutation, repeat until converged
			// Pick the best one
			// Stop once:
			// 1. The csm is less than 1e-4
			// 2. The difference between the old and new csm is less than 1%
			// 3. The max number of iterations has been reached
			while ((fabs(dist) > 1e-4) && (fabs(old - dist)/fabs(old) > 0.01) && (iterNum < maxIters)) {
				old = dist;
				estimatePerm(m, temp, tempDir, type);
				runSinglePerm(m, outAtoms, temp, &dist, tempDir, dMin, type);					
				if (dist < best) {
					best = dist;
					memcpy(bestPerm, temp, sizeof(int) * m->_size);
					memcpy(bestDir, tempDir, sizeof(double) * 3);
				}	
				iterNum++;			
				//printf("Old csm: %6.4f New csm %6.4f\n", old, dist);
			};

			// Keep the best solution so far...
			if (best < *csm) { 
				*csm = best;
				memcpy(perm, bestPerm , sizeof(int) * m->_size);
				memcpy(dir, bestDir, sizeof(double) * 3);
			}		
			printf("Attempt #%d: best csm is %4.2f after %d iterations\n", (i+1), best, iterNum);			
		}
		for (i = 0; i < n_dirs; i++) {
			free(dirs[i]);
		}
		free(dirs);
		free(bestPerm);	

		// run once more to get the out atoms right !
		runSinglePerm(m, outAtoms, perm, csm, dir, dMin, type);					
	}
	
	free(temp);	
	
}

/*
 * Find an initial guess for the approximated symmetry direction, 
 * which in the case of mirror symmetry is vec perpendicular to the mirror plane, 
 * and in other cases, the symmetry axis.
 * This is done using least-mean-squares, which provides 3 guesses, times 3 if we try to remove outliers
 */
void findSymmetryDirection(Molecule *m, double  ***dirs, int *n_dirs, OperationType type) {
	int* groupSizes = (int*)malloc(sizeof(int) * m->_groupNum);
	double** groupAverages = (double**)malloc(sizeof(double*) * m->_groupNum);
	int *outliers = (int*)malloc(sizeof(int) * m->_groupNum);
	int i,j;
	double **testDir;
	double median;
	double zero[3] = {0.0,0.0,0.0};
	std::vector<double *> results;
	int useOrthogonal = TRUE;

	testDir = (double**)malloc(sizeof(double*) * 3);
	for (i = 0; i < 3; i++) {
		testDir[i] = (double*)malloc(sizeof(double) * 3);
	}	
		
	for (i = 0; i < m->_groupNum; i++) { 
		groupSizes[i] = 0;
		groupAverages[i] = (double*)malloc(sizeof(double) * 3);
		groupAverages[i][0] = groupAverages[i][1] = groupAverages[i][2] = 0.0;
	}
	
	for (i = 0; i < m->_size; i++) { 
		groupSizes[m->_similar[i] - 1]++;
		groupAverages[m->_similar[i] - 1][0] += m->_pos[i][0];
		groupAverages[m->_similar[i] - 1][1] += m->_pos[i][1];
		groupAverages[m->_similar[i] - 1][2] += m->_pos[i][2];
	}

	for (i = 0; i < m->_groupNum; i++) { 					
		outliers[i] = FALSE;
		groupAverages[i][0] /= 	groupSizes[i];	
		groupAverages[i][1] /= 	groupSizes[i];	
		groupAverages[i][2] /= 	groupSizes[i];	
	}

	if (type == CS) { 
		// For CS 			
		// Assuming that all orbit center-of-masses are 
		// found on the plane of reflection - find it			
		// for some reason - we get a few options
		planeFit(groupAverages, m->_groupNum, testDir,outliers);		
	} else {
		// For CN and SN with N > 2			
		// Assuming that all orbit center-of-masses are found on the axis of symmetry - find it			
		// for some reason - we get a few options
		lineFit(groupAverages, m->_groupNum, testDir, outliers);		
	}	

	// if there are not enough groups for reliable outlier detection - don't invoke it.
	if (detectOutliers && m->_groupNum >= MIN_GROUPS_FOR_OUTLIERS) {			
		for (j = 0; j < 3; j++) { 
			double **tempDir;

			tempDir = (double**)malloc(sizeof(double*) * 3);
			for (i = 0; i < 3; i++) {
				tempDir[i] = (double*)malloc(sizeof(double) * 3);
			}	

			// 1. Find the distance of each point from the line / plane
			// 2. Find the median m
			// 3. for each distance di, if (di / m > A || m / di > A - remove as outlier)
			// 4. recompute line / plane
			double* dists = (double *)malloc(m->_groupNum * sizeof(double));				
			for (i = 0; i < m->_groupNum; i++) {
				if (type == CS) { 
					dists[i] = fabs(testDir[j][0] * groupAverages[i][0] + 
						testDir[j][1] * groupAverages[i][1] + 
						testDir[j][2] * groupAverages[i][2]);							
				} else { 
					dists[i] = computeDistanceFromLine(groupAverages[i],zero,testDir[j]);
				}

			}
			median = findMedian(dists, m->_groupNum);
			for (i = 0; i < m->_groupNum; i++) {
				if (dists[i] / median > A || dists[i] / median > A) {					
					outliers[i] = true; 
				} 
			}
			if (type == CS) {	
				planeFit(groupAverages, m->_groupNum, tempDir, outliers);		
			} else {
				lineFit(groupAverages, m->_groupNum, tempDir, outliers);		
			}
			for (i = 0; i < 3; i++) { 	
				results.push_back(tempDir[i]);
			}

			free(dists);		
			free(tempDir);		
		}
		for (i = 0; i < 3; i++) {
			free(testDir[i]);
		}			

	} else {					
		// just copy...
		for (i = 0; i < 3; i++) {	
			results.push_back(testDir[i]);
		}
		
	}

	if (useOrthogonal) {	
		// Go over all dirs and add two orthogonal axes
		std::vector<double*> temp = results;
		for (i = 0; i < temp.size(); i++) { 
			double* dir = temp[i];
			double* newDir = (double *)malloc(sizeof(double) * 3);
			double* newDir2 = (double *)malloc(sizeof(double) * 3);
			double norm, scal;
			if (fabs(dir[0]) < MINDOUBLE) { 
				newDir[0] = 1.0; newDir[1] = 0.0; newDir[2] = 0.0;
				newDir2[0] = 0.0; newDir2[1] = -dir[2]; newDir2[2] = dir[1];
			} else if (fabs(dir[1]) < MINDOUBLE) {
				newDir[0] = -dir[2]; newDir[1] = 0.0; newDir[2] = dir[0];
				newDir2[0] = newDir2[2] = 0.0; newDir[1] = 1.0;	
			} else {
				newDir[0] = -dir[1]; newDir[1] = dir[0]; newDir[2] = 0.0;
				newDir2[0] = 0.0; newDir2[1] = -dir[2]; newDir2[2] = dir[1];
			}
			// normalize first
			norm = sqrt(newDir[0] * newDir[0] + newDir[1] * newDir[1] + newDir[2] * newDir[2]);
			newDir[0] /= norm; newDir[1] /= norm; newDir[2] /= norm;
			// remove projection of first from second (both already orthogonal to dir
			scal = newDir[0] * newDir2[0] + newDir[1] * newDir2[1] + newDir[2] * newDir2[2];
			newDir2[0] -= scal * newDir[0]; newDir2[1] -= scal * newDir[1]; newDir2[2] -= scal * newDir[2];
			// normalize second
			norm = sqrt(newDir2[0] * newDir2[0] + newDir2[1] * newDir2[1] + newDir2[2] * newDir2[2]);
			newDir2[0] /= norm; newDir2[1] /= norm; newDir2[2] /= norm;

			results.push_back(newDir); 
			results.push_back(newDir2);

		}
	}



	// initialize results array
	*n_dirs = results.size();
	(*dirs) = (double**)malloc(sizeof(double*) *(*n_dirs));
	for (i = 0; i < *n_dirs; i++) {				
		(*dirs)[i] = results[i];
	}

	for (i = 0; i < m->_groupNum; i++) { 		
		free(groupAverages[i]);
	}
	free(outliers);
	free(groupSizes);
	free(groupAverages);
	free(testDir);
}

void estimatePerm(Molecule* m, int *perm, double *dir, OperationType type) {
	int isImproper = (type != CN) ? TRUE : FALSE;
	int isZeroAngle = (type == CS) ? TRUE : FALSE;
	int maxGroupSize = getMaxGroupSize(m);
	int *group = (int*)malloc(sizeof(int) * maxGroupSize);
	int *used = (int*)malloc(sizeof(int) * m->_size);
	int i, j, k, l;
	double rotaionMatrix[3][3];
	double tmpMatrix[3][3] = {{0.0, -dir[2], dir[1]}, {dir[2], 0.0, -dir[0]}, {-dir[1], dir[0], 0.0}};
	double angle;
	double **rotated = (double**)malloc(sizeof(double*) * m->_size);
	struct distRecord * distances = (struct distRecord *)malloc(sizeof(struct distRecord) * maxGroupSize * maxGroupSize);
	int factor = (isImproper ? (-1) : 1);	
	int tableSize = 0;
	int orbitDone, orbitSize, orbitStart;
	int left;
	angle = isZeroAngle ? 0.0 : (2 * PI / opOrder);
		
	for (j = 0; j < m->_size; j++) {
		rotated[j] = (double *)malloc(sizeof(double) * 3);		
		rotated[j][0] = rotated[j][1] = rotated[j][2] = 0;
		perm[j] = -1;		
	}		

	// Prepare the rotation matrix
	for (j = 0; j < 3; j++) {
		for (k = 0; k < 3; k++) {
			rotaionMatrix[j][k] = 
				((j == k) ? cos(angle) : 0) + 
				(factor - cos(angle)) * dir[j] * dir[k] + 
				sin(angle) * tmpMatrix[j][k];
		}		
	}

	// Run the operation on the current point set
	for (j = 0; j < m->_size; j++) {
		for (k = 0; k < 3; k++) {
			for (l = 0; l < 3; l++) {
				rotated[j][k] += rotaionMatrix[k][l] * m->_pos[j][l];				
			}
		}
		//printf("%d (%4.2f, %4.2f, %4.2f) -> (%4.2f, %4.2f, %4.2f)\n", j,
		//	m->_pos[j][0], m->_pos[j][1], m->_pos[j][2], 
		//	rotated[j][0], rotated[j][1], rotated[j][2]);
	}

	// run over the groups
	for (i = 0; i < m->_groupNum; i++) { 
		// Get the group
		int groupSize = getGroup(m, i + 1,group);		
		for (j = 0; j < m->_size; j++) { 
			used[j] = 0;
		}
		
		// compute the distance matrix
		for (j = 0; j <	groupSize; j++) { 
			for (k = 0; k < groupSize; k++) { 	
					int index = j * groupSize + k;
					distances[index].row = group[j];
					distances[index].col = group[k];
					distances[index].distance = 0;
					for (l = 0; l < 3; l++ ) { 
						distances[index].distance += 
							(m->_pos[group[j]][l] - rotated[group[k]][l]) * 
							(m->_pos[group[j]][l] - rotated[group[k]][l]);
					}		
					distances[index].distance = sqrt(distances[index].distance);
			}		
		}
	
		tableSize = groupSize * groupSize;

		// Sort the distances			
		qsort(distances, tableSize, sizeof(struct distRecord), distComp);
	
		/*	
		printf("Working on group: ");
		for (j = 0; j < groupSize; j++) { 
			printf("%d ", group[j]);
		}
		printf("\n");
		*/
		

		left = groupSize;			
		// Go over the sorted group, and set the permutation
		for (j = 0; j < tableSize && left > 0; j++) { 			
			int enoughForFullOrbit = left >= opOrder;			
			int row = distances[j].row;
			int col = distances[j].col;				
	
			// If we have used this item already - skip it.
			if (perm[row] != -1) continue;							
			//printf("%d %d -- %f\n", row, col,distances[j].distance);
			
			// If we do not have enought to full groups, set all remaining items to themselves
			if (left == 1 || (type == CN && !enoughForFullOrbit)) { 							
				for (k = 0; k < groupSize; k++) { 
					if (used[group[k]] == 0) { 
						perm[group[k]] = group[k];						
						//printf("set %d<->%d\n", group[k], group[k]);
					}				
				}				
				break;
			} 			
			if (opOrder == 2) {			
				// Special treatment - only size 1 and two orbits are allowed 
				// If both elements are not yet set, use the element.			
				if (perm[row] == -1 && perm[col] == -1)  { 
					perm[row] = col;
					perm[col] = row;					
					//printf("set %d<->%d\n", row, col);
					left -= (row == col) ? 1 : 2;
				}
			} else {
				// we now want to complete an orbit. 
				if (perm[row] == -1 && used[col] == 0) {					
					perm[row] = col;
					used[col] = 1;				
					//printf("set %d->%d\n", row, col);
					left--;						
				} else {
					continue;
				}

				// if this is an orbit of size one...
				if (row == col) continue;

				// If there is no more room for full orbit, must be SN and size two orbit
				if (type == SN && !enoughForFullOrbit) { 
					perm[col] = row;
					used[row] = 1;					
					//printf("set %d->%d\n", col, row);
					left--;					
					continue;
				}

				// Run until an orbit is complete
				orbitDone = FALSE;
				orbitStart = row;				
				orbitSize = 1;	
				
				while (!orbitDone) {									
					if (orbitSize == opOrder - 1) { 
						//printf("Closing orbit...\n");
						//Close the group - we've reached the end							
						row = col;
						col = orbitStart;
						orbitDone = TRUE;
					} else {								
						// Search for the next orbit element
						for (k = j + 1; k < tableSize; k++) { 
							if (distances[k].row == col && used[distances[k].col] == 0 && 			
								distances[k].col != distances[k].row) {
								if (orbitStart == distances[k].col) { 
									if (type == SN && orbitSize ==1) { 
										// we have now closed an orbit of size 2			
										orbitDone = TRUE;
									} else {	
										continue;
									}
								}
								row = distances[k].row;
								col = distances[k].col;
								orbitSize++;								
								break;		
							}
						}								
					}									
					perm[row] = col;
					used[col] = 1;
					//printf("set %d->%d\n", row, col);	
					left--;					
				}
			}		
		} 	
						
	}	

	// verify that the orbits are correct?
	
	for (j = 0; j < m->_size; j++) { 		
		free(rotated[j]);			
	}	

	free(rotated);	
	free(distances);
	free(group);
}

void lineFit(double **points, int nPoints, double **dirs, int *outliers) {
	// taken from http://www.mapleprimes.com/forum/linear-regression-in-3d
	double A[3] = {0,0,0};	
	double matrix[3][3] = { {0,0,0}, {0,0,0}, {0,0,0}};

	double **copyMat = dmatrix(1,3,1,3);	
	double *diag = dvector(1,3);
	double *secdiag = dvector(1,3);	
	int realNum = 0;

	double norm;
	int i,j,k;	

	// Compute A
	for (i = 0 ; i < nPoints; i++) { 
		if (!outliers[i]) {
			A[0] += points[i][0];
			A[1] += points[i][1];
			A[2] += points[i][2];
			realNum ++;
		}
	}
	A[0] /= realNum;
	A[1] /= realNum;
	A[2] /= realNum;
	
	// Compute matrix
	for (i = 0; i < nPoints; i++) { 
		if (!outliers[i]) {
			for (j = 0; j < 3; j++) { 
				for (k = 0; k < 3; k++) { 				
					matrix[j][k] += (points[i][j] - A[j]) * (points[i][k] - A[k]);
				}
			}
		}
	}

	// perhaps actual copying is needed
	for (i = 0; i < 3; i++) {		
		for (j = 0; j < 3; j++) {
			copyMat[i + 1][j + 1] = matrix[i][j];
		}
	}		

	// compute the matrix's eigenvalues and eigenvectors.
	// In the end - diag is eigenvals
	// cols of matrix are eigen vecs
	tred2(copyMat, 3, diag, secdiag);
	tqli(diag, secdiag, 3, copyMat);

	// We just try the three lines?...	
	for (i = 0; i < 3; i++) {		
		for (j = 0; j < 3; j++) {
			dirs[i][j] = copyMat[j + 1][i + 1];
		}
	}

	for (i = 0; i < 3; i++) { 
		norm = sqrt(dirs[i][0] * dirs[i][0] + dirs[i][1] * dirs[i][1] + dirs[i][2] * dirs[i][2]); 
		dirs[i][0] /= norm;
		dirs[i][1] /= norm;
		dirs[i][2] /= norm;

	}	

	free_dmatrix(copyMat, 1, 3, 1, 3);	
	free_dvector(diag, 1, 3);
	free_dvector(secdiag, 1, 3);	
}
void planeFit(double **points, int nPoints, double **dirs, int* outliers) {
	// taken from http://www.mapleprimes.com/forum/linear-regression-in-3d
	double A[3] = {0,0,0};	
	double matrix[3][3] = { {0,0,0}, {0,0,0}, {0,0,0}};

	double **copyMat = dmatrix(1,3,1,3);	
	double *diag = dvector(1,3);
	double *secdiag = dvector(1,3);	

	double norm;
	int i,j,k;
	int realNum = 0;
	
	// Compute A
	for (i = 0 ; i < nPoints; i++) { 
		if (!outliers[i]) {
			realNum++;
			A[0] += points[i][0];
			A[1] += points[i][1];
			A[2] += points[i][2];
		}
	}
	A[0] /= realNum;
	A[1] /= realNum;
	A[2] /= realNum;

	// Compute matrix
	for (i = 0; i < nPoints; i++) { 
		if (outliers[i]) continue;
		for (j = 0; j < 3; j++) { 
			for (k = 0; k < 3; k++) { 
				matrix[j][k] += (points[i][j] - A[j]) * (points[i][k] - A[k]);
			}
		}
	}

	// perhaps actual copying is needed
	for (i = 0; i < 3; i++) {		
		for (j = 0; j < 3; j++) {
			copyMat[i + 1][j + 1] = matrix[i][j];
		}
	}		

	// compute the matrix's eigenvalues and eigenvectors.
	// In the end - diag is eigenvals
	// cols of matrix are eigen vecs
	tred2(copyMat, 3, diag, secdiag);
	tqli(diag, secdiag, 3, copyMat);

	// We just try the three planes...	
	for (i = 0; i < 3; i++) {		
		for (j = 0; j < 3; j++) {
			dirs[i][j] = copyMat[j + 1][i + 1];
		}		
	}	
	for (i = 0; i < 3; i++) { 
		norm = sqrt(dirs[i][0] * dirs[i][0] + dirs[i][1] * dirs[i][1] + dirs[i][2] * dirs[i][2]); 
		dirs[i][0] /= norm;
		dirs[i][1] /= norm;
		dirs[i][2] /= norm;

	}	

	free_dmatrix(copyMat, 1, 3, 1, 3);	
	free_dvector(diag, 1, 3);
	free_dvector(secdiag, 1, 3);	
}


/*
* prints the Molecule position, outcome position, csm, dMin and directional cosines to output file
*/
void printOutput(Molecule* m, double** outAtoms, double csm, double *dir, double dMin, FILE *out, double* localCSM){

	int i,j;
	printf("%s: %.6lf\n",opName,fabs(csm));
	fprintf(out, "%s: %.4lf\n",opName,fabs(csm));
	fprintf(out, "SCALING FACTOR: %7lf\n", dMin);

	fprintf(out, "\n INITIAL STRUCTURE COORDINATES\n%i\n",m->_size);

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

	fprintf(out, "\n RESULTING STRUCTURE COORDINATES\n%i\n",m->_size);

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

	if (printNorm) {
		printf( "NORMALIZATION FACTOR: %7lf\n", m->_norm);
		printf( "SCALING FACTOR OF SYMMETRIC STRUCTURE: %7lf\n", dMin);
		printf( "DIRECTIONAL COSINES: %lf %lf %lf\n", dir[0], dir[1], dir[2]);
		printf( "NUMBER OF EQUIVALENCE GROUPS: %d\n", m->_groupNum);
	}

	if (printLocal) {
		double sum = 0;
		fprintf(out,"\nLocal CSM: \n");	
		for (i = 0; i < m->_size; i++) {
			sum += localCSM[i];
			fprintf(out,"%s %7lf\n", m->_symbol[i], localCSM[i]);
		}
		fprintf(out,"\nsum: %7lf\n", sum);
	}


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

	if (printNorm) {
		printf( "NORMALIZATION FACTOR: %7lf\n", m->_norm);
		printf( "SCALING FACTOR OF SYMMETRIC STRUCTURE: %7lf\n", dMin);
		printf( "DIRECTIONAL COSINES: %lf %lf %lf\n", dir[0], dir[1], dir[2]);
		printf( "NUMBER OF EQUIVALENCE GROUPS: %d\n", m->_groupNum);
	}

}

/*
* prints in PDB format the Molecule position, outcome position, csm, dMin and directional cosines to output file
*/
void printOutputFormat(Molecule* m, OBMol& mol, double** outAtoms, double csm, double *dir, double dMin, FILE *out, char *fname, double* localCSM) {

	fprintf(out, "%s: %.4lf\n",opName,fabs(csm));
	fprintf(out, "SCALING FACTOR: %7lf\n", dMin);

	// TODO - should we print the centered molecule, or the original one (and, accordingly, the symmetric struct)

	fprintf(out, "\n INITIAL STRUCTURE COORDINATES\n");

	updateCoordinates(mol, m->_pos);

	writeMolecule(mol, format, out, fname);

	updateCoordinates(mol, outAtoms);

	fprintf(out, "\n RESULTING STRUCTURE COORDINATES\n");

	writeMolecule(mol, format, out, fname);

	// print results to screen

	fprintf(out, "\n DIRECTIONAL COSINES:\n\n");
	fprintf(out, "%lf %lf %lf\n", dir[0], dir[1], dir[2]);

	if (writeOpenu)
		printf("SV* %.4lf *SV\n",fabs(csm));
	else
		printf( "%s: %.4lf\n",opName,fabs(csm));

	if (printNorm) {
		printf( "NORMALIZATION FACTOR: %7lf\n", m->_norm);
		printf( "SCALING FACTOR OF SYMMETRIC STRUCTURE: %7lf\n", dMin);
		printf( "DIRECTIONAL COSINES: %lf %lf %lf\n", dir[0], dir[1], dir[2]);
		printf( "NUMBER OF EQUIVALENCE GROUPS: %d\n", m->_groupNum);
	}

	if (printLocal) {
		double sum = 0;
		int i;
		fprintf(out,"\nLocal CSM: \n");	
		for (i = 0; i < m->_size; i++) {
			sum += localCSM[i];
			fprintf(out,"%s %7lf\n", m->_symbol[i], localCSM[i]);
		}
		fprintf(out,"\nsum: %7lf\n", sum);
	}
}

/*------------------------------------------------------------------------*/
/*---------------------------   MARK       -------------------------------*/

/*------    Pinsky -------------------*/


/*
 * Normalizes the position of atoms of the molecule
 * returns one [TRUE] if successful, zero[FALSE] otherwise
 */
int normalizeMolecule2(Molecule *m){

	double tmp,x_avg, y_avg, z_avg,norm;
	int i;

	x_avg = y_avg = z_avg = 0.0;

	for(i=0; i< m->_size; i++){
		x_avg += m->_pos[i][0];
		y_avg += m->_pos[i][1];
		z_avg += m->_pos[i][2];
	}
	x_avg /= (double)(m->_size);
	y_avg /= (double)(m->_size);
	z_avg /= (double)(m->_size);

	norm = 0.0;
	for(i=0; i< m->_size; i++){
		tmp = SQR(m->_pos[i][0]-x_avg) +
		      SQR(m->_pos[i][1]-y_avg) +
		      SQR(m->_pos[i][2]-z_avg);
		norm += tmp;
	}

	norm = sqrt(norm / (double)m->_size);
	

	if(norm < MINDOOUBLE)
		return FALSE;

	for(i=0; i< m->_size; i++){
		m->_pos[i][0] = ((m->_pos[i][0] - x_avg) / norm);
		m->_pos[i][1] = ((m->_pos[i][1] - y_avg) / norm);
		m->_pos[i][2] = ((m->_pos[i][2] - z_avg) / norm);
	}

	return TRUE;
}



#ifdef CNV

double calc_cnv(Molecule* m, double **rotAtoms, double *dir, int nSize,
              int **perm_matrix)
{

	int i, ii, j, k, l;
	
	double angle, B1, B2, det, en_value_min, sum, sym, cn_sum;

	double rotaionMatrix[3][3];
	
	double a_matrix[3][3]  = {{0.0, 0.0, 0.0}, 
	                          {0.0, 0.0, 0.0}, 
				  {0.0, 0.0, 0.0}};
	
	double tmpMatrix[3][3] = {{0.0, -dir[2], dir[1]}, 
	                          {dir[2], 0.0, -dir[0]}, 
				  {-dir[1], dir[0], 0.0}};
				  
	double EMatrix[3][3]   = {{1.0, 0.0, 0.0}, 
	                          {0.0, 1.0, 0.0}, 
				  {0.0, 0.0, 1.0}};
				  
	double **copyMat = dmatrix(1,3,1,3);
	double *diag = dvector(1,3);
	double *secdiag = dvector(1,3);			  
				  
	double temp[3],scalar[3];
	

        /*    Calculation of A-matrix    */
						  				  
		  		  	
	for (i=0; i < nSize; i++) {
	
	    angle=i*PI/nSize;
	    
	    
	    for (j = 0; j < 3; j++)
	        for (k = 0; k < 3; k++)
		    rotaionMatrix[j][k] = cos(angle)*EMatrix[j][k]+ 
				      (1.0-cos(angle))*dir[j]*dir[k]+ 
				      sin(angle)*tmpMatrix[j][k];
				      
	    for (j = 0; j < m->_size; j++)
	        for (k = 0; k < 3; k++)
	            rotAtoms[j][k] = 0.0;
		
	    for (j = 0; j < m->_size; j++)
	        for (k = 0; k < 3; k++)
		    for (l = 0; l < 3; l++)
		        rotAtoms[j][k] += rotaionMatrix[k][l] * m->_pos[j][l];
	     	     	
	     for (j=0; j < m->_size; j++) {		
	     	ii = perm_matrix[nSize+i][j];
		     
	        for (k = 0; k < 3; k++)
		    for (l = 0; l < 3; l++)	     
	               a_matrix[k][l] += rotAtoms[j][k]*rotAtoms[ii][l]+	     		
		                         rotAtoms[ii][k]*rotAtoms[j][l];		 	
			
	    }
			
	}
		
	
	/*     Calculation of Cn transformation      */
	

	cn_sum = 0.0;
	for (i=0; i < nSize; i++) {
	
	    angle=i*2.0*PI/nSize;
	    

	    for (j = 0; j < 3; j++)
	        for (k = 0; k < 3; k++)
		    rotaionMatrix[j][k] = cos(angle)*EMatrix[j][k]+ 
				      (1.0-cos(angle))*dir[j]*dir[k]+ 
				      sin(angle)*tmpMatrix[j][k];		
				      
	    for (j = 0; j < m->_size; j++)
	        for (k = 0; k < 3; k++)
	            rotAtoms[j][k] = 0.0;	
	
	    for (j = 0; j < m->_size; j++)
	        for (k = 0; k < 3; k++)
		    for (l = 0; l < 3; l++)
		        rotAtoms[j][k] += rotaionMatrix[k][l] * m->_pos[j][l];
			
	     for (j=0; j < m->_size; j++) {		
	     	 ii = perm_matrix[i][j];
	         for (k = 0; k < 3; k++)
		       cn_sum += m->_pos[j][k]*rotAtoms[ii][k];
		       
	     }				
        }	
	
	
		
	/* Calculation of eigenvalues and eigenvectors of A_matrix  */	
	/* Calculation of scalar multiplications of eigenvectors with axes   */	
	

	for (i = 0; i < 3; i++)
	    for (j = 0; j < 3; j++)
		copyMat[i + 1][j + 1] = a_matrix[i][j];		

	tred2(copyMat, 3, diag, secdiag);
	
	tqli(diag, secdiag, 3, copyMat);
	
	for (i = 0; i < 3; i++) {
	    scalar[i] = 0.0;
	    for (j = 0; j < 3; j++)
	       scalar[i] += dir[j]*copyMat[j+1][i+1];	       
	    temp[i] = scalar[i] * scalar[i];
	}
	
	
	/* Solution of quadratic equation  */
	
	
	B1 = temp[0]*(diag[2]+diag[3])+
	     temp[1]*(diag[1]+diag[3])+
	     temp[2]*(diag[1]+diag[2]);
	    
	B2 = temp[0]*(diag[2]*diag[3])+
	     temp[1]*(diag[1]*diag[3])+
	     temp[2]*(diag[1]*diag[2]);
	     
	det = 0.25*B1*B1-B2;
	
	
	if (det < 0.0) return(MAXDOUBLE);	
	
	en_value_min = 0.5*B1-sqrt(det);

		
	/* Calculation of Cnv measure   */


        sum = 0.0;	
	for (i=0; i < nSize; i++)
	    for (j = 0; j < m->_size; j++)	    
	        for (k = 0; k < 3; k++)  {
		   ii   = perm_matrix[nSize+i][j]; 
	           sum += m->_pos[j][k]*m->_pos[ii][k];
		   		   	   
	        }
				
												
/*        sym = 1.0+(en_value_min-sum)/(m->_size*nSize);*/

        sym = 1.0+(en_value_min-sum-cn_sum)/(2.0*m->_size*nSize);

	free_dmatrix(copyMat, 1, 3, 1, 3);
	free_dvector(diag, 1, 3);
	free_dvector(secdiag, 1, 3);
	
	return(sym);		    		
}





#else // DN
double calc_dn(Molecule* m, double **rotAtoms, double *dir, int nSize,
              int **perm_matrix)
{

	int i, ii, j, k, l;
	
	double angle, B1, B2, det, en_value_min, en_value_max, sum, cn_sum, sym;

	double rotaionMatrix[3][3];
	
	double a_matrix[3][3]  = {{0.0, 0.0, 0.0}, 
	                          {0.0, 0.0, 0.0}, 
				  {0.0, 0.0, 0.0}};
	
	double tmpMatrix[3][3] = {{0.0, -dir[2], dir[1]}, 
	                          {dir[2], 0.0, -dir[0]}, 
				  {-dir[1], dir[0], 0.0}};
				  
	double EMatrix[3][3]   = {{1.0, 0.0, 0.0}, 
	                          {0.0, 1.0, 0.0}, 
				  {0.0, 0.0, 1.0}};
				  
	double **copyMat = dmatrix(1,3,1,3);
	double *diag = dvector(1,3);
	double *secdiag = dvector(1,3);			  
				  
	double temp[3],scalar[3];
	

        /*    Calculation of A-matrix    */
						  				  
		  		  	
	for (i=0; i < nSize; i++) {
	
	    angle=i*PI/nSize;
	    
	    
	    for (j = 0; j < 3; j++)
	        for (k = 0; k < 3; k++)
		    rotaionMatrix[j][k] = cos(angle)*EMatrix[j][k]+ 
				      (1.0-cos(angle))*dir[j]*dir[k]+ 
				      sin(angle)*tmpMatrix[j][k];
				      
	    for (j = 0; j < m->_size; j++)
	        for (k = 0; k < 3; k++)
	            rotAtoms[j][k] = 0.0;
		
	    for (j = 0; j < m->_size; j++)
	        for (k = 0; k < 3; k++)
		    for (l = 0; l < 3; l++)
		        rotAtoms[j][k] += rotaionMatrix[k][l] * m->_pos[j][l];
	     	     	
	     for (j=0; j < m->_size; j++) {		
	     	ii = perm_matrix[nSize+i][j];
		     
	        for (k = 0; k < 3; k++)
		    for (l = 0; l < 3; l++)	     
	               a_matrix[k][l] += rotAtoms[j][k]*rotAtoms[ii][l]+	     		
		                         rotAtoms[ii][k]*rotAtoms[j][l];		 	
			
	    }
			
	}
	
	
	
	/*     Calculation of Cn transformation      */
	

	cn_sum = 0.0;
	for (i=0; i < nSize; i++) {
	
	    angle=i*2.0*PI/nSize;
	    

	    for (j = 0; j < 3; j++)
	        for (k = 0; k < 3; k++)
		    rotaionMatrix[j][k] = cos(angle)*EMatrix[j][k]+ 
				      (1.0-cos(angle))*dir[j]*dir[k]+ 
				      sin(angle)*tmpMatrix[j][k];		
				      
	    for (j = 0; j < m->_size; j++)
	        for (k = 0; k < 3; k++)
	            rotAtoms[j][k] = 0.0;	
	
	    for (j = 0; j < m->_size; j++)
	        for (k = 0; k < 3; k++)
		    for (l = 0; l < 3; l++)
		        rotAtoms[j][k] += rotaionMatrix[k][l] * m->_pos[j][l];
			
	     for (j=0; j < m->_size; j++) {		
	     	 ii = perm_matrix[i][j];
	         for (k = 0; k < 3; k++)
		       cn_sum += m->_pos[j][k]*rotAtoms[ii][k];
		       
	     }				
        }	
		
	
		
	/* Calculation of eigenvalues and eigenvectors of A_matrix  */	
	/* Calculation of scalar multiplications of eigenvectors with axes   */	
	

	for (i = 0; i < 3; i++)
	    for (j = 0; j < 3; j++)
		copyMat[i + 1][j + 1] = a_matrix[i][j];		

	tred2(copyMat, 3, diag, secdiag);
	
	tqli(diag, secdiag, 3, copyMat);
	
	for (i = 0; i < 3; i++) {
	    scalar[i] = 0.0;
	    for (j = 0; j < 3; j++)
	       scalar[i] += dir[j]*copyMat[j+1][i+1];	       
	    temp[i] = scalar[i] * scalar[i];
	}
	
	
	/* Solution of quadratic equation  */
	
	
	B1 = temp[0]*(diag[2]+diag[3])+
	     temp[1]*(diag[1]+diag[3])+
	     temp[2]*(diag[1]+diag[2]);
	    
	B2 = temp[0]*(diag[2]*diag[3])+
	     temp[1]*(diag[1]*diag[3])+
	     temp[2]*(diag[1]*diag[2]);
	     
	det = 0.25*B1*B1-B2;
	
	
	if (det < 0.0) return(MAXDOUBLE);	
	
/*	en_value_min = 0.5*B1-sqrt(det);*/
	en_value_max = 0.5*B1+sqrt(det);

		
	/* Calculation of Dn measure   */


        sum = 0.0;	
	for (i=0; i < nSize; i++)
	    for (j = 0; j < m->_size; j++)	    
	        for (k = 0; k < 3; k++)  {
		   ii   = perm_matrix[nSize+i][j]; 
	           sum += m->_pos[j][k]*m->_pos[ii][k];		   
		   
	        }
		
												
/*        sym = 1.0+(en_value_min-sum)/(m->_size*nSize);
	sym = 1.0+(sum-en_value_max)/(m->_size*nSize); */
	
        sym = 1.0+(sum-cn_sum-en_value_max)/(2.0*m->_size*nSize);	


	free_dmatrix(copyMat, 1, 3, 1, 3);
	free_dvector(diag, 1, 3);
	free_dvector(secdiag, 1, 3);
	
	return(sym);		    		
}
#endif
