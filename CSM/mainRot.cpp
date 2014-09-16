/*
 * Author: Amir Zait
 *
 * This is the main program. 
 * It includes the algorithm itself, as well as dealing with the program's IO
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> //for strcmp,strlen
#include "groupPermuter.h"
#include <openbabel/mol.h>
#include "babelAdapter.h"
#include <vector>
#include "Molecule.h"

#include "logging.h"
#include <boost/format.hpp>
#include <iomanip>
#include <sstream>

using namespace std;

extern "C" {
#include "mainhelpers.h"
#include "nrutil.h"
}

#define CSMFORMAT "CSM"
#define MAXDOUBLE  100000000.0
#define MINDOUBLE  1e-8
#define GROUPSIZE_LIMIT 15
#define GROUPSIZE_FACTOR 1.32e11
#define APPROX_RUN_PER_SEC 8e4
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
void findBestPermUsingDir(Molecule* m, double** outAtoms, int* optimalPerm, double* csm, double* dir, double* dMin, OperationType type);
void findSymmetryDirection(Molecule *m, double  ***dirs, int *n_dirs, OperationType type);
void estimatePerm(Molecule* m, int *perm, double *dir, OperationType type);
void readPerm(FILE* permfile, int* perm, int size);
void readDir(FILE* dirFile, double* dir);
void lineFit(double **points, int nPoints, double **dirs, int* outlies);
void planeFit(double **points, int nPoints, double **dirs, int *outliers);
void initIndexArrays(Molecule* m, int* posToIdx, int* idxToPos);
double createSymmetricStructure(Molecule* m, double **outAtom, int *perm, double *dir, OperationType type, double dMin);
double computeLocalCSM(Molecule* m, double *localCSM, int *perm, double *dir, OperationType type);

// global options
bool ignoreHy = false;
bool removeHy = false;
bool ignoreSym = false;
bool useFormat = false;
bool writeOpenu = false;
OperationType type;
int opOrder;
bool useperm = false;
bool useDir = false;
bool findPerm = false;
bool useMass = false;
bool limitRun = true;
char *format = NULL;
bool babelBond = false;
bool timeOnly = false;
int sn_max = 4;
//int anneal = false;
bool detectOutliers = false;
double A = 2;
bool babelTest = false;
bool printNorm = false;
bool printLocal = false;
bool keepCenter = false;
bool debug = false;

// file pointers
FILE* inFile = NULL;
FILE* outFile = NULL;
FILE* permfile = NULL;
FILE* dirfile = NULL;
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
	printf("-debug			- Enable debug printings");
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
		for (i = 1; i <= m->groupNum(); i++) {			
			int groupSize = m->getGroupSize(i);
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
		for (i = 1; i <= m->groupNum(); i++) {
			int groupSize = m->getGroupSize(i);
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
 * returns one [true] if successful, zero[false] otherwise
 */
void normalize(double **coords, Molecule *m){

	double tmp,x_avg, y_avg, z_avg,norm;
	int i;

	x_avg = y_avg = z_avg = 0.0;

	if (!keepCenter) {
		double mass_sum = 0;
		for(i=0; i< m->size(); i++){
			x_avg += coords[i][0] * m->mass(i);
			y_avg += coords[i][1] * m->mass(i);
			z_avg += coords[i][2] * m->mass(i);
			mass_sum += m->mass(i);
		}
		x_avg /= (double)(mass_sum);
		y_avg /= (double)(mass_sum);
		z_avg /= (double)(mass_sum);
	}

	norm = 0.0;
	for(i=0; i< m->size(); i++){
		tmp = SQR(coords[i][0]-x_avg) +
		      SQR(coords[i][1]-y_avg) +
		      SQR(coords[i][2]-z_avg);
		norm += tmp;
	}
	// normalize to 1 and not molecule size
	//norm = sqrt(norm / (double)m->size());
	norm = sqrt(norm);


	for(i=0; i< m->size(); i++){
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
			LOG(fatal) << "Only Even values of n are allowed";
			exit(1);
		}
		sprintf(opName, "S%d SYMMETRY",opOrder);	
	}
	// try to open infile for reading
	inFileName = argv[2];
	if ((inFile = fopen(inFileName, "rt")) == NULL){
		if (writeOpenu) 
		{
			printf("ERR* Failed to open data file %s *ERR\n", inFileName);
		}
		LOG(fatal) << "Failed to open data file " << inFileName;
		exit(1);
	}


	// try to open outfile for writing
	outFileName = argv[3];
	if ((outFile = fopen(outFileName, "w")) == NULL){
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
	for ( i=4;  i< argc ;  i++ ){
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
		} else if (nextIsDirFile) {
			char* dirfilename = argv[i];
			if ((dirfile = fopen(dirfilename, "rt")) == NULL){
				if (writeOpenu) {
					printf("ERR* Failed to open dir file %s for reading *ERR\n", dirfilename);
				}
				LOG(fatal) << "Failed to open dir file " << dirfilename << " for reading";
				exit(1);
			}
			nextIsDirFile = false;
		} else if (nextIsMaxSn) { 
			sn_max = atoi(argv[i]);
			nextIsMaxSn = false;
		} else if (strcmp(argv[i],"-sn_max" ) == 0) {
			if (type != CH) { 
				LOG(fatal) << "Option -sn_max only applies to chirality";
				exit(1);
			}
			nextIsMaxSn = true;	
		} else if (strcmp(argv[i],"-ignoreHy" ) == 0 )
			ignoreHy = true;
		else if (strcmp(argv[i],"-removeHy" ) == 0 )
			removeHy = true;

		else if (strcmp(argv[i],"-ignoreSym" ) == 0 )
			ignoreSym = true;

	    else if (strncmp(argv[i],"-format", 7 ) == 0 ) {
		useFormat = true;
		format = strdup(argv[i] + 7);
		} else if (strcmp(argv[i],"-writeOpenu" ) == 0 ) {
			writeOpenu = true;	
		} else if (strcmp(argv[i], "-nolimit") == 0) { 
			limitRun = false;  
		} else if (strcmp(argv[i], "-useperm") == 0) {
			useperm = true;
			nextIsPermFile = true;
		} else if (strcmp(argv[i], "-usedir") == 0) {
			useDir = true;
			nextIsDirFile = true;
		} else if (strcmp(argv[i], "-babelbond") == 0) {
			babelBond = true;
		} else if (strcmp(argv[i], "-useMass") == 0) { 
			useMass = true;					
		} else if (strcmp(argv[i], "-timeonly") == 0) {
			timeOnly = true;
		} else if  (strcmp(argv[i], "-printNorm") == 0) {
			printNorm = true;
		} else if (strcmp(argv[i], "-help") == 0) {
			usage(argv[0]);
			exit(0);
		} else if (strcmp(argv[i], "-findperm") == 0) { 
			findPerm = true;
		} else if (strcmp(argv[i], "-detectOutliers") == 0) {
			detectOutliers = true;
		} else if (strcmp(argv[i], "-approx") == 0) {
			detectOutliers = true;
			findPerm = true;
		} else if (strcmp(argv[i], "-babelTest") == 0) { 
			babelTest = true;
		} else if (strcmp(argv[i], "-printlocal") == 0) { 
			printLocal = true;
		} else if (strcmp(argv[i], "-keepCenter") == 0) { 
			keepCenter = true;
		} else if (strcmp(argv[i], "-debug") == 0) {
			debug = true;
		}
	}
	if (writeOpenu) {
		useFormat = true;
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
	double *localCSM = NULL;		 

	init_logging();

	LOG(info) << "CSM starting up";

	// init options
	parseInput(argc,argv);
	set_debug_logging(debug);

	if ((findPerm && useperm) || (findPerm && useDir) || (useDir && useperm)) { 
		LOG(fatal) << "-findperm, -useperm and -usedir are mutually exclusive";
		exit(1);
	} 

	// try to read molecule from infile
	Molecule* m;
	OBMol mol; 
	OperationType chMinType = CS;
	int chMinOrder = 2;
	
	if (useFormat) {
		// If a specific format is used, read molecule using that format
		if (strcasecmp(format, CSMFORMAT) == 0) {
			m = Molecule::create(inFile,stdout,ignoreSym && !useperm);
			if (m==NULL) exit(1);
			if (useMass)
			{
				m->fillAtomicMasses();
			}
		} else {
			mol = readMolecule (inFileName, format, babelBond);
			m = Molecule::createFromOBMol(mol, ignoreSym && !useperm, useMass);	
		}
   	} else {
		format = getExtension(inFileName);

		// if the extension is CSM - use csm
		if (strcasecmp(format, CSMFORMAT) == 0) {
			m = Molecule::create(inFile,stdout,ignoreSym && !useperm);
			if (m==NULL) exit(1);
			if (useMass)
			{
				m->fillAtomicMasses();
			}
		} else {
			
			mol = readMolecule (inFileName, NULL, babelBond);
			m = Molecule::createFromOBMol(mol, ignoreSym && !useperm, useMass);						
		}
   	}

	if (babelTest) // Mol is ok - return 0
		return 0;

	if (!m){
		if (writeOpenu) {
			printf("ERR* Failed to read molecule from data file *ERR\n");
		}
		LOG(fatal) << "Failed to read molecule from data file";
		exit(1);
	}

	// strip unwanted atoms if needbe
	if ((ignoreHy || removeHy) && !useperm){
		char* removeList[] = {"H"," H"};
		Molecule* n = NULL;
		if (ignoreHy)
			n = m->stripAtoms(removeList,2,false);
		else //removeHy 
			n = m->stripAtoms(removeList,2,true);		
	
		mol.DeleteHydrogens();
		
	
		if (!n){
			if (writeOpenu) {
				printf("ERR* Failed while trying to strip unwanted atoms *ERR\n");
			}
			LOG(fatal) << "Failed while trying to strip unwanted atoms";
			exit(1);
		}
		delete m;
		m = n;
		// continue as per usual
	}

	if (!findPerm) {
		if (!useperm && !useDir) {
			double time = 1.0*totalNumPermutations(m) / 3600 / APPROX_RUN_PER_SEC;
			if (time != time) {
				// time is NaN
				time = MAXDOUBLE;
			}
			LOG(info) << "Going to enumerate over " << totalNumPermutations(m) << " permutations";
			LOG(info) << "Entire run should take approx. " << std::setprecision(2) << std::fixed << time << " hours on a 2.0Ghz Computer";
		} else {
			LOG(info) << "Using 1 permutation";
			LOG(info) << "Run should be instantaneous";
		} 
		if (timeOnly) { return 0; };
	}

	// allocate memory for outAtoms
	outAtoms = (double **)malloc(m->size() * sizeof(double*));
	for (i=0;i<m->size();i++)
		outAtoms[i] = (double *)malloc(3 * sizeof(double));
       
	perm = (int *)malloc(sizeof(int) * m->size());

	//normalize Molecule
	if (!m->normalizeMolecule(keepCenter)){
		if (writeOpenu) {
			printf("ERR* Failed to normalize atom positions: dimension of set of points = zero *ERR\n");
		}
		LOG(fatal) << "Failed to normalize atom positions: dimension of set of points = zero";
		exit(1);
	}
 
	if (useDir) readDir(dirfile,dir);

	if (useperm) {	
		if (type == CH) {
			LOG(fatal) << "Chirality can't be given a permutation, run the specific csm operation instead";
			exit(1);
		}	
		readPerm(permfile,perm, m->size());
		runSinglePerm(m, outAtoms, perm, &csm, dir, &dMin, type);
	} else {
		if (type != CH) { 
			// perform operation
			if (useDir) {
				findBestPermUsingDir(m, outAtoms, perm, &csm, dir, &dMin, type);				
			} else if (findPerm) { 
				findBestPerm(m, outAtoms, perm, &csm, dir, &dMin, type);				
			} else {
				csmOperation(m, outAtoms, perm, &csm, dir, &dMin, type);			
			}
		} else {
			// chirality support 
			double chCsm, chdMin;
			double **chOutAtoms;
			int *chPerm;			
			double chDir[3] = {0.0, 0.0, 0.0};						

			chOutAtoms = (double **)malloc(m->size() * sizeof(double*));
			for (i=0;i<m->size();i++)
				chOutAtoms[i] = (double *)malloc(3 * sizeof(double));
       
			chPerm = (int *)malloc(sizeof(int) * m->size());
			chMinType = CS;						
			opOrder = 2;		
			if (useDir) {
				findBestPermUsingDir(m, outAtoms, perm, &csm, dir, &dMin, CS);				
			} else if (findPerm) { 
				findBestPerm(m, outAtoms, perm, &csm, dir, &dMin, CS);				
			} else { 
				csmOperation(m, outAtoms, perm, &csm, dir, &dMin, CS);			
			}

			if (csm > MINDOUBLE) {							
				for (i = 2; i <= sn_max; i+=2) {				
					opOrder = i;
					if (useDir) {
						findBestPermUsingDir(m, chOutAtoms, chPerm, &chCsm, chDir, &chdMin, SN);
					} else if (findPerm) { 
						findBestPerm(m, chOutAtoms, chPerm, &chCsm, chDir, &chdMin, SN);
					} else {
						csmOperation(m, chOutAtoms, chPerm, &chCsm, chDir, &chdMin, SN);
					}				
					if (chCsm < csm) {
						int j;
						chMinType = SN;
						chMinOrder = opOrder;
						csm = chCsm;
						dMin = chdMin;
						memcpy(dir, chDir, sizeof(double)* 3);
						memcpy(perm, chPerm, sizeof(int) * m->size());	
						for (j = 0; j < m->size(); j++) { 	
							memcpy(outAtoms[j], chOutAtoms[j], sizeof(double)*3);
						}
					}
					if (csm < MINDOUBLE) break;
				}
			}
			// housekeeping
			for (i=0;i<m->size();i++){
				free(chOutAtoms[i]);
			}
			free(chOutAtoms);	
			free(chPerm);
		}		
	}

	if (printLocal) {	
		localCSM = (double *)malloc(sizeof(double) * m->size());
		if (type == CH) opOrder = chMinOrder;
		computeLocalCSM(m,localCSM, perm, dir,  type != CH ? type : chMinType);
	}

	normalize(outAtoms, m);

	// De-normalize
	for (i = 0; i < m->size(); i++) { 
		m->pos()[i][0] *= m->norm();
		m->pos()[i][1] *= m->norm();
		m->pos()[i][2] *= m->norm();
		outAtoms[i][0] *= m->norm();
		outAtoms[i][1] *= m->norm();
		outAtoms[i][2] *= m->norm();
	}	

   	if (useFormat) {
		// If a specific format is used, read molecule using that format
		if (strcasecmp(format, CSMFORMAT) == 0) {
			printOutput(m, outAtoms, csm, dir, dMin, outFile, localCSM);
		} else {
			printOutputFormat(m, mol, outAtoms, csm, dir, dMin, outFile, outFileName, localCSM);
		}
	} else {
		// if the extension is CSM - use csm
		if (strcasecmp(getExtension(inFileName), CSMFORMAT) == 0) {
			printOutput(m, outAtoms, csm, dir, dMin, outFile, localCSM);
		} else {
			printOutputFormat(m, mol, outAtoms, csm, dir, dMin, outFile, outFileName, localCSM);
		}
	}

	if (type == CH) {
		if (chMinType == CS) { 		
			fprintf(outFile, "\n MINIMUM CHIRALITY WAS FOUND IN CS\n\n");
		} else { 
			fprintf(outFile, "\n MINIMUM CHIRALITY WAS FOUND IN S%d\n\n", chMinOrder);
		}
	}	

	fprintf(outFile, "\n PERMUTATION:\n\n");
	for (i = 0; i < m->size(); i++) {
		fprintf(outFile, "%d ", perm[i] + 1);
	}
	fprintf(outFile,"\n");

	// housekeeping
	for (i=0;i<m->size();i++){
		free(outAtoms[i]);
	}
	free(outAtoms);
	delete m;
	free(perm);

	if (printLocal) free(localCSM);

	fclose(inFile);
	fclose(outFile);

	if (permfile != NULL)
		fclose(permfile);

	if (dirfile != NULL) 
		fclose(dirfile);

	return 0;
}

void readDir(FILE* dirfile, double* dir) { 
	fscanf(dirfile, "%lf%lf%lf", &dir[0], &dir[1], &dir[2]);
}

/**
 * Read a permutation from a file - simply a list separated by spaces.
 */
void readPerm(FILE* permfile, int* perm, int size) {
	int *used = (int *)malloc(sizeof(int) * size);
	int i = 0;
	for (i = 0; i < size; i++) {
		used[i] = false;
	}
	for (i = 0; i < size; i++) {
		int cur = -1;
		int res = fscanf(permfile,"%d", &cur);
		if (res != 1 || cur < 1 || cur > size || used[cur - 1]) {
			if (writeOpenu) {
				printf("ERR* Invalid permutation *ERR\n");
			}
			LOG(fatal) << "Invalid permutation";
			free(used);
			fclose(permfile);
			exit(1);
		}
		used[cur - 1] = true;
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
	for ( j=1;  j<= m->groupNum() ;  j++ ){
		for ( i=0;  i< m->size() ;  i++ ){
			if (m->similar(i) == j){
				idxToPos[counter] = i;
				counter++;
			}
		}
	}

	// build posToIdx
	for ( i=0;  i< m->size() ;  i++ ){
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
	double **copyMat = dmatrix(1, 3, 1, 3);
	double *copyVec = dvector(1, 3);
	double *diag = dvector(1, 3);
	double *secdiag = dvector(1, 3);
	double *temp = dvector(1, 3);

	double matrix[3][3] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };
	double vec[3] = { 0.0, 0.0, 0.0 };
	double scalar[3];
	double maxval, scl, angle;
	double rtr[6] { 0, 0, 0, 0, 0, 0 }, rti[6] { 0, 0, 0, 0, 0, 0 };
	double coeffs[7];
	int i,j;
	int *curPerm = (int *)malloc(sizeof(int) * m->size());
	double csm, dists;

	int isImproper = (type != CN) ? true : false;
	int isZeroAngle = (type == CS) ? true : false;

	LOG(debug) << "calcRefPlane called";
	stringstream permstrm;
	for (int i = 0; i < m->size(); i++)
	{
		if (i > 0)
			permstrm << ", ";
		permstrm << perm[i];
	}
	LOG(debug) << "Permutation is " << permstrm.str();
	LOG(debug) << "Direction is " << setprecision(2) << fixed << dir[0] << " " << dir[1] << " " << dir[2];

	// initialize identity permutation
	for (i = 0; i < m->size(); i++) {
		curPerm[i] = i;
	}

	// compute matrices according to current perm and its powers (the identity does not contribute anyway)
	for (i = 1; i < opOrder; i++) {			
		angle = isZeroAngle ? 0.0 : (2 * PI * i / opOrder);
		// The i'th power of the permutation
		for (j = 0; j < m->size(); j++) {
			curPerm[j] = perm[curPerm[j]];
		}
		if (isImproper && ((i % 2) == 1)) {
			computeMatrix(m->pos(), curPerm, m->size(), &matrix, -1-cos(angle));	
		} else {
			computeMatrix(m->pos(), curPerm, m->size(), &matrix, 1-cos(angle));	
		}	  
		computeVector(m->pos(), curPerm, m->size(), &vec, sin(angle));
	}

	LOG(debug).unsetf(ios_base::fixed);
	LOG(debug) << "Computed matrix is:" << setprecision(4);
	LOG(debug) << matrix[0][0] << " " << matrix[0][1] << " " << matrix[0][2];
	LOG(debug) << matrix[1][0] << " " << matrix[1][1] << " " << matrix[1][2];
	LOG(debug) << matrix[2][0] << " " << matrix[2][1] << " " << matrix[2][2];

	// perhaps actual copying is needed
	for (i = 0; i < 3; i++) {
		copyVec[i + 1] = vec[i];
		for (j = 0; j < 3; j++) {
			copyMat[i + 1][j + 1] = matrix[i][j];
		}
	}		

	LOG(debug) << "Copied matrix is:";
	LOG(debug) << copyMat[1][1] << " " << copyMat[1][2] << " " << copyMat[1][3];
	LOG(debug) << copyMat[2][1] << " " << copyMat[2][2] << " " << copyMat[2][3];
	LOG(debug) << copyMat[3][1] << " " << copyMat[3][2] << " " << copyMat[3][3];


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

	LOG(debug) <<  "copyVec: " << copyVec[1] << " " << copyVec[2] << " " << copyVec[3];

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
	LOG(debug) << "Coefficients: " << coeffs[0] << ", " << coeffs[1] << ", " << coeffs[2] << ", " << coeffs[3] << ", " << coeffs[4] << ", " << coeffs[5] << ", " << coeffs[6];
	rpoly(coeffs, 6, rtr, rti);

	LOG(debug) << "rtr: " << rtr[0] << " " << rtr[1] << " " << rtr[2] << " " << rtr[3] << " " << rtr[4] << " " << rtr[5];
	LOG(debug) << "rti: " << rti[0] << " " << rti[1] << " " << rti[2] << " " << rti[3] << " " << rti[4] << " " << rti[5];

	maxval = -MAXDOUBLE;
	for (i = 0; i < 6; i++) {								
		if (maxval < rtr[i] && (fabs(rti[i]) < ZERO_IM_PART_MAX)) maxval = rtr[i];		
	}	

	LOG(debug) << setprecision(6) << fixed << "diag: " << diag[1] << " " << diag[2] << " " << diag[3];
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
	for (i = 0; i < m->size(); i++) {
		curPerm[i] = i;
	}

	// compute CSM	
	csm = 0.0;		
	for (i = 0; i < opOrder; i++) {		
		// This can be more efficient - save results of matrix computation
		if (i != 0) {
			angle = isZeroAngle ? 0.0 : ((2 * PI * i) / opOrder);
			dists = 0.0;
			for (j = 0; j < m->size(); j++) {
				// i'th power of permutation
				curPerm[j] = perm[curPerm[j]];	
			}
			for (j = 0; j < m->size(); j++) {
				dists += (m->pos()[j][0] * m->pos()[curPerm[j]][0] + 
					m->pos()[j][1] * m->pos()[curPerm[j]][1] + 
					m->pos()[j][2] * m->pos()[curPerm[j]][2]); 
			}
			csm += cos(angle) * dists;

		} else {
			csm += 1.0;	
		}
	}	

	
	LOG(debug) << setprecision(6) << fixed << "csm=" << csm << " maxval=" << maxval << " scl=" << scl;
	LOG(debug) << setprecision(6) << fixed << "dir: " << dir[0] << " " << dir[1] << " " << dir[2];
	csm += (maxval -  scl) / 2;
	csm = fabs(100 * (1.0 - csm / opOrder));
	free(curPerm);

	LOG(debug) << setprecision(6) << fixed << "dir - csm: " << dir[0] << " " << dir[1] << " " << dir[2] << " - " << csm;
	
	free_dmatrix(copyMat, 1, 3, 1, 3);
	free_dvector(copyVec, 1, 3);
	free_dvector(diag, 1, 3);
	free_dvector(secdiag, 1, 3);
	free_dvector(temp,1,3);
	return csm;
}

double createSymmetricStructure(Molecule* m, double **outAtoms, int *perm, double *dir, OperationType type, double dMin) {
	int isImproper = (type != CN) ? true : false;
	int isZeroAngle = (type == CS) ? true : false;
	int i, j, k, l;
	int *curPerm = (int *)malloc(sizeof(int) * m->size());
	double rotaionMatrix[3][3];
	double tmpMatrix[3][3] = {{0.0, -dir[2], dir[1]}, {dir[2], 0.0, -dir[0]}, {-dir[1], dir[0], 0.0}};
	double angle;
	double res = 0.0;

	for (i = 0; i < m->size(); i++) {
		// initialize with identity operation
		curPerm[i] = i;
		for (j = 0; j < 3; j++) {
			outAtoms[i][j] = m->pos()[i][j];
		}
	}

	for (i = 1; i < opOrder; i++) {
		angle = isZeroAngle ? 0.0 : (2 * PI * i / opOrder);
		int factor = ((isImproper && (i % 2) == 1) ? (-1) : 1);
		for (j = 0; j < m->size(); j++) {
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

		for (j = 0; j < m->size(); j++) {
			for (k = 0; k < 3; k++) {
				for (l = 0; l < 3; l++) {
					outAtoms[j][k] += rotaionMatrix[k][l] * m->pos()[curPerm[j]][l];
				}
			}			
		}

	}

	for (j = 0; j < m->size(); j++) {
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
	int isImproper = (type != CN) ? true : false;
	int isZeroAngle = (type == CS) ? true : false;
	int i, j, k, l;
	int *curPerm = (int *)malloc(sizeof(int) * m->size());
	double rotaionMatrix[3][3];
	double tmpMatrix[3][3] = {{0.0, -dir[2], dir[1]}, {dir[2], 0.0, -dir[0]}, {-dir[1], dir[0], 0.0}};
	double angle;
	double res = 0.0;

	double tempCSM = 0.0;


	for (i = 0; i < m->size(); i++) {
		// initialize with identity operation
		curPerm[i] = i;
		localCSM[i] = 0;
	}

	for (i = 1; i < opOrder; i++) {
		angle = isZeroAngle ? 0.0 : (2 * PI * i / opOrder);
		int factor = ((isImproper && (i % 2) == 1) ? (-1) : 1);
		for (j = 0; j < m->size(); j++) {
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

		for (j = 0; j < m->size(); j++) {
			double rotated[3] = {0.0,0.0,0.0};
			for (k = 0; k < 3; k++) {
				for (l = 0; l < 3; l++) {
					rotated[k] += rotaionMatrix[k][l] * m->pos()[curPerm[j]][l];
				}
			}			
			for (k = 0; k < 3; k++) {
				tempCSM += SQR(rotated[k] - m->pos()[j][k]);
				localCSM[j] += SQR(rotated[k] - m->pos()[j][k]);
			}			
		}

	}

	for (j = 0; j < m->size(); j++) {
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
	double curDir[3] = { 0, 0, 0 };
	int *idxToPos, *posToIdx;
	int * groupSizes;
	GroupPermuter* gp;
	int addGroupsOfTwo;

	// These are the 
	int *realPerm = (int *)malloc(sizeof(int) * m->size());

	// allocate memory for index arrays arrays
	idxToPos = (int*)malloc(m->size() * sizeof(int));
	posToIdx = (int*)malloc(m->size() * sizeof(int));

	// and group sizes arrays
	groupSizes = (int*)malloc(m->groupNum() * sizeof(int*));	

	for(i=0; i<m->size(); i++)    	
		// init index arrays
		initIndexArrays(m,posToIdx,idxToPos);

	// init csm, curCsm
	*csm = curCsm = MAXDOUBLE;

	// get group sizes
	for(i=1; i<= m->groupNum() ; i++){
		groupSizes[i-1] = m->getGroupSize(i);
	}

	// create permuter
	if (type == SN && opOrder > 2) {
		addGroupsOfTwo = 1;
	} else {
		addGroupsOfTwo = 0;
	}
	gp = new GroupPermuter(m->groupNum(),groupSizes,m->size(),opOrder, addGroupsOfTwo);
	if (!gp){
		if (writeOpenu) {
			printf("ERR* Failed to create groupPermuter *ERR\n");
		}
		LOG(fatal) << "Failed to create groupPermuter";
		exit(1);
	};

	// calculate csm for each valid permutation & remember minimal (in optimalAntimer)
	while ( gp->next() ) {

		for (i = 0; i < m->size(); i++) {
			realPerm[i] = idxToPos[gp->elementAt(posToIdx[i])];
		}			
		curCsm = calcRefPlane(m, realPerm, curDir, type);		

		// check, if it's a minimal csm, update maxGroupCsm and optimalAntimer
		if(curCsm < *csm) {
			*csm = curCsm;
			dir[0] = curDir[0]; dir[1] = curDir[1]; dir[2] = curDir[2];
			
			for (i = 0; i < m->size(); i++) {
				optimalPerm[i] = realPerm[i];
			}
		}				
	}

	// failed to find value for any permutation
	
	if (*csm == MAXDOUBLE){
		if (writeOpenu) {
			printf("ERR* Failed to calculate a csm value for %s *ERR\n",opName);
		}
		LOG(fatal) << "Failed to calculate a csm value for " << opName;
		exit(1);
	}

	// which is DMIN?
	*dMin = (1.0 - (*csm / 100 * opOrder / (opOrder - 1)));
	createSymmetricStructure(m, outAtoms, optimalPerm, dir, type, *dMin);
		

	// housekeeping
	free(groupSizes);
	free(idxToPos);
	free(posToIdx);
	delete gp;
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

void findBestPermUsingDir(Molecule* m, double** outAtoms, int* perm, double* csm, double* dir, double* dMin, OperationType type) {
	estimatePerm(m, perm, dir, type);
	runSinglePerm(m, outAtoms, perm, csm, dir, dMin, type);	
}

/**
 * Finds an approximate permutation which can be used in the analytical computation.
 */
void findBestPerm(Molecule* m, double** outAtoms, int* perm, double* csm, double* dir, double* dMin, OperationType type) {
	LOG(debug) << "findBestPerm called";

	int *temp = (int*)malloc(sizeof(int) * m->size());
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
		int *bestPerm = (int*)malloc(sizeof(int) * m->size());	
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
					memcpy(bestPerm, temp, sizeof(int) * m->size());
					memcpy(bestDir, tempDir, sizeof(double) * 3);
				}	
				iterNum++;			
				LOG(debug) << "Old csm: " << setprecision(4) << fixed << old << ", new csm " << dist;
			};

			// Keep the best solution so far...
			if (best < *csm) { 
				*csm = best;
				memcpy(perm, bestPerm , sizeof(int) * m->size());
				memcpy(dir, bestDir, sizeof(double) * 3);
			}		
			LOG(info) << "Attempt #" << i + 1 << ": best csm is " << setprecision(2) << fixed << best << " after " << iterNum << " iterations";
		}
		for (i = 0; i < n_dirs; i++) {
			free(dirs[i]);
		}
		free(dirs);
		free(bestPerm);	

		// run once more to get the out atoms right !
		runSinglePerm(m, outAtoms, perm, csm, dir, dMin, type);					
	}
		
#if 0 
	if (anneal) {
			
		// The initial probability is set to be 2 / N 
		// for a change of 0.0001
		double initialTemp = 0.0001 * log(1.0 * m->size()/4.0);
		// printf("IT = %4.2f\n", initialTemp);
		double alpha = 0.9;
		int stepsPerPart = 2000;
		double t;			
		double dist;
		double old = *csm;
		int factor = 1;
		double flipProb = 0.001;

		int **groups = (int**)malloc(sizeof(int *) * m->groupNum());
		int *groupSizes = (int*) malloc(sizeof(int) * m->groupNum());
		for (i = 0; i < m->groupNum(); i++ ) {
			groupSizes[i] = getGroupSize(m, i + 1); 
			groups[i] = (int*)malloc(sizeof(int) * groupSizes[i]);
			getGroup(m, i + 1, groups[i]);
		}

		srand( time (NULL) );
#ifndef _WIN32
		srand48( time (NULL) );			
#endif
		// try to anneal the best permutation
		memcpy(temp, perm, sizeof(int) * m->size());
		for (t = initialTemp; t >= 0.0001*initialTemp; t *= alpha) { 
			for (i = 0; i < stepsPerPart /* * m->size() */; i++) { 
				// select a node at random
				int first = rand() % m->size();

				// select another from its similarity group
				int second = groups
					[m->similar(first) - 1]
					[rand() % groupSizes[m->similar(first) - 1]];
				
				int temp1 = first, temp2 = second, firstSize = 1, secondSize = 1;
				if (first == second) continue;
					
				// find the orbit size of the two nodes, as well as the atoms preceeding them	
				while (temp[temp1] != first) {	
					firstSize++;
					temp1 = temp[temp1];
				}	
				while (temp[temp2] != second) {	
					secondSize++;
					temp2 = temp[temp2];
				}	
				// Now, finally, operate
				if (firstSize == 1 && secondSize == 1) { 		
					// both are single-orbits - so do nothing
				} else if (firstSize == 1) {
					// only first is a single orbit
					temp[temp2] = first;
					temp[first] = temp[second];
					temp[second] = second;
				} else if (secondSize == 1) { 
					// only second is a single orbit
					temp[temp1] = second;
					temp[second] = temp[first];
					temp[first] = first;
				} else {
					// both are not single orbits
					int t = temp[first];	
					temp[first] = temp[second];
					temp[second] = t;	
					t = temp[temp1];
					temp[temp1] = temp[temp2];
					temp[temp2] = t;
				}
				
				runSinglePerm(m, outAtoms, temp, &dist, dir, dMin, type);	
				if (factor == 1) { 
					if (drand48() < flipProb) { 
						factor = -1;
					}
				} else { 	
					// spend one tenth of the time in reverted state
					if (drand48() < (flipProb * 100)) { 	
						factor = 1;
					}
				}
				// printf("Tried to change from %4.2f to %4.2f: ", factor * old, factor * dist);
				if (dist < old || drand48() < exp(factor * (old - dist) / t)) {	
					// if this improves - 
					
					/*if (dist < old) { 
						printf("Better\n");	
					} else {
						printf("Kept\n");
					} */
					if (dist < *csm) { 
						*csm = dist;
						printf("Changed to %4.2f\n", *csm);
						memcpy(perm, temp , sizeof(int) * m->size());
					}
					old = dist;
				} else {
						// undo			
					//printf("Reverted\n");								
					if (firstSize == 1 && secondSize == 1) { 		
						// both are single-orbits - so do nothing
					} else if (firstSize == 1) {
						// only first is a single orbit
						temp[temp2] = second;
						temp[second] = temp[first];
						temp[first] = first;
					} else if (secondSize == 1) { 
						// only second is a single orbit
						temp[temp1] = first;
						temp[first] = temp[second];
						temp[second] = second;
					} else {
						// both are not single orbits
						int t = temp[first];	
						temp[first] = temp[second];
						temp[second] = t;	
						t = temp[temp1];
						temp[temp1] = temp[temp2];
						temp[temp2] = t;
					}
					
				} 
			}	
			
			// printf("At T = %4.2f - csm %4.2f\n", t, *csm);
		}

		for (i = 0; i < m->groupNum(); i++) {
			free(groups[i]);
		}
		free(groups);
		free(groupSizes);
	}					
#endif

	free(temp);	
	
}

/*
 * Find an initial guess for the approximated symmetry direction, 
 * which in the case of mirror symmetry is vec perpendicular to the mirror plane, 
 * and in other cases, the symmetry axis.
 * This is done using least-mean-squares, which provides 3 guesses, times 3 if we try to remove outliers
 */
void findSymmetryDirection(Molecule *m, double  ***dirs, int *n_dirs, OperationType type) {
	int* groupSizes = (int*)malloc(sizeof(int) * m->groupNum());
	double** groupAverages = (double**)malloc(sizeof(double*) * m->groupNum());
	int *outliers = (int*)malloc(sizeof(int) * m->groupNum());
	int i,j;
	double **testDir;
	double median;
	double zero[3] = {0.0,0.0,0.0};
	std::vector<double *> results;
	int useOrthogonal = true;

	testDir = (double**)malloc(sizeof(double*) * 3);
	for (i = 0; i < 3; i++) {
		testDir[i] = (double*)malloc(sizeof(double) * 3);
	}	
		
	for (i = 0; i < m->groupNum(); i++) { 
		groupSizes[i] = 0;
		groupAverages[i] = (double*)malloc(sizeof(double) * 3);
		groupAverages[i][0] = groupAverages[i][1] = groupAverages[i][2] = 0.0;
	}
	
	for (i = 0; i < m->size(); i++) { 
		groupSizes[m->similar(i) - 1]++;
		groupAverages[m->similar(i) - 1][0] += m->pos()[i][0];
		groupAverages[m->similar(i) - 1][1] += m->pos()[i][1];
		groupAverages[m->similar(i) - 1][2] += m->pos()[i][2];
	}

	for (i = 0; i < m->groupNum(); i++) { 					
		outliers[i] = false;
		groupAverages[i][0] /= 	groupSizes[i];	
		groupAverages[i][1] /= 	groupSizes[i];	
		groupAverages[i][2] /= 	groupSizes[i];	
	}

	if (type == CS) { 
		// For CS 			
		// Assuming that all orbit center-of-masses are 
		// found on the plane of reflection - find it			
		// for some reason - we get a few options
		planeFit(groupAverages, m->groupNum(), testDir,outliers);		
	} else {
		// For CN and SN with N > 2			
		// Assuming that all orbit center-of-masses are found on the axis of symmetry - find it			
		// for some reason - we get a few options
		lineFit(groupAverages, m->groupNum(), testDir, outliers);		
	}	

	// if there are not enough groups for reliable outlier detection - don't invoke it.
	if (detectOutliers && m->groupNum() >= MIN_GROUPS_FOR_OUTLIERS) {			
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
			double* dists = (double *)malloc(m->groupNum() * sizeof(double));				
			for (i = 0; i < m->groupNum(); i++) {
				if (type == CS) { 
					dists[i] = fabs(testDir[j][0] * groupAverages[i][0] + 
						testDir[j][1] * groupAverages[i][1] + 
						testDir[j][2] * groupAverages[i][2]);							
				} else { 
					dists[i] = computeDistanceFromLine(groupAverages[i],zero,testDir[j]);
				}

			}
			median = findMedian(dists, m->groupNum());
			for (i = 0; i < m->groupNum(); i++) {
				if (dists[i] / median > A || dists[i] / median > A) {					
					outliers[i] = true; 
				} 
			}
			if (type == CS) {	
				planeFit(groupAverages, m->groupNum(), tempDir, outliers);		
			} else {
				lineFit(groupAverages, m->groupNum(), tempDir, outliers);		
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

	for (i = 0; i < m->groupNum(); i++) { 		
		free(groupAverages[i]);
	}
	free(outliers);
	free(groupSizes);
	free(groupAverages);
	free(testDir);
}

void estimatePerm(Molecule* m, int *perm, double *dir, OperationType type) {
	int isImproper = (type != CN) ? true : false;
	int isZeroAngle = (type == CS) ? true : false;
	int maxGroupSize = m->getMaxGroupSize();
	int *group = (int*)malloc(sizeof(int) * maxGroupSize);
	int *used = (int*)malloc(sizeof(int) * m->size());
	int i, j, k, l;
	double rotaionMatrix[3][3];
	double tmpMatrix[3][3] = {{0.0, -dir[2], dir[1]}, {dir[2], 0.0, -dir[0]}, {-dir[1], dir[0], 0.0}};
	double angle;
	double **rotated = (double**)malloc(sizeof(double*) * m->size());
	struct distRecord * distances = (struct distRecord *)malloc(sizeof(struct distRecord) * maxGroupSize * maxGroupSize);
	int factor = (isImproper ? (-1) : 1);	
	int tableSize = 0;
	int orbitDone, orbitSize, orbitStart;
	int left;
	angle = isZeroAngle ? 0.0 : (2 * PI / opOrder);

	LOG(debug) << "estimatePerm called";

	for (j = 0; j < m->size(); j++) {
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
	for (j = 0; j < m->size(); j++) {
		for (k = 0; k < 3; k++) {
			for (l = 0; l < 3; l++) {
				rotated[j][k] += rotaionMatrix[k][l] * m->pos()[j][l];				
			}
		}
		LOG(debug) << boost::format("%d (%4.2f, %4.2f, %4.2f) -> (%4.2f, %4.2f, %4.2f)") %
			j % 
			m->pos()[j][0] % m->pos()[j][1] % m->pos()[j][2] % 
			rotated[j][0] % rotated[j][1] % rotated[j][2];
	}

	// run over the groups
	for (i = 0; i < m->groupNum(); i++) { 
		// Get the group
		int groupSize = m->getGroup(i + 1,group);		
		for (j = 0; j < m->size(); j++) { 
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
							(m->pos()[group[j]][l] - rotated[group[k]][l]) * 
							(m->pos()[group[j]][l] - rotated[group[k]][l]);
					}		
					distances[index].distance = sqrt(distances[index].distance);
			}		
		}
	
		tableSize = groupSize * groupSize;

		// Sort the distances			
		qsort(distances, tableSize, sizeof(struct distRecord), distComp);
	
		stringstream groupstrm;
		for (j = 0; j < groupSize; j++)
			groupstrm << group[j] << " ";
		LOG(debug) << "Working on group " << groupstrm.str();
		

		left = groupSize;			
		// Go over the sorted group, and set the permutation
		for (j = 0; j < tableSize && left > 0; j++) { 			
			int enoughForFullOrbit = left >= opOrder;			
			int row = distances[j].row;
			int col = distances[j].col;				
	
			// If we have used this item already - skip it.
			if (perm[row] != -1) 
				continue;							
			LOG(debug) << row << " " << col << " " << distances[j].distance;
			
			// If we do not have enought to full groups, set all remaining items to themselves
			if (left == 1 || (type == CN && !enoughForFullOrbit)) { 							
				for (k = 0; k < groupSize; k++) { 
					if (used[group[k]] == 0) { 
						perm[group[k]] = group[k];						
						LOG(debug) << "set " << group[k] << "<->" << group[k];
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
					LOG(debug) << "set " << row << "<->" << col;
					left -= (row == col) ? 1 : 2;
				}
			} else {
				// we now want to complete an orbit. 
				if (perm[row] == -1 && used[col] == 0) {					
					perm[row] = col;
					used[col] = 1;				
					LOG(debug) << "set " << row << "<->" << col;
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
					LOG(debug) << "set " << row << "->" << col;
					left--;
					continue;
				}

				// Run until an orbit is complete
				orbitDone = false;
				orbitStart = row;				
				orbitSize = 1;	
				
				while (!orbitDone) {									
					if (orbitSize == opOrder - 1) { 
						LOG(debug) << "Closing orbit";
						row = col;
						col = orbitStart;
						orbitDone = true;
					} else {								
						// Search for the next orbit element
						for (k = j + 1; k < tableSize; k++) { 
							if (distances[k].row == col && used[distances[k].col] == 0 && 			
								distances[k].col != distances[k].row) {
								if (orbitStart == distances[k].col) { 
									if (type == SN && orbitSize ==1) { 
										// we have now closed an orbit of size 2			
										orbitDone = true;
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
					LOG(debug) << "set " << row << "->" << col;
					left--;
				}
			}		
		} 	
						
	}	

	// verify that the orbits are correct?
	
	for (j = 0; j < m->size(); j++) { 		
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

	fprintf(out, "\n INITIAL STRUCTURE COORDINATES\n%i\n",m->size());

	for(i=0; i<m->size(); i++){
		fprintf(out, "%3s%10lf %10lf %10lf\n",
			m->symbol(i), m->pos()[i][0], m->pos()[i][1], m->pos()[i][2]);
	}

	for (i = 0; i < m->size(); i++) {
		fprintf(out, "%d ", i + 1);
		for ( j = 0; j < m->valency(i); j++ ) {
			fprintf(out, "%d ", m->adjacent(i,j) + 1);
		}
		fprintf(out,"\n");
	}

	fprintf(out, "\n RESULTING STRUCTURE COORDINATES\n%i\n",m->size());

	for(i=0; i<m->size(); i++){
		fprintf(out, "%3s%10lf %10lf %10lf\n",
			m->symbol(i), outAtoms[i][0], outAtoms[i][1], outAtoms[i][2]);
	}

	for (i = 0; i < m->size(); i++) {
		fprintf(out, "%d ", i + 1);
		for ( j = 0; j < m->valency(i); j++ ) {
			fprintf(out, "%d ", m->adjacent(i,j) + 1);
		}
		fprintf(out,"\n");
	}

	fprintf(out, "\n DIRECTIONAL COSINES:\n\n");
	fprintf(out, "%lf %lf %lf\n", dir[0], dir[1], dir[2]);

	if (printNorm) {
		printf( "NORMALIZATION FACTOR: %7lf\n", m->norm());
		printf( "SCALING FACTOR OF SYMMETRIC STRUCTURE: %7lf\n", dMin);
		printf( "DIRECTIONAL COSINES: %lf %lf %lf\n", dir[0], dir[1], dir[2]);
		printf( "NUMBER OF EQUIVALENCE GROUPS: %d\n", m->groupNum());
	}

	if (printLocal) {
		double sum = 0;
		fprintf(out,"\nLocal CSM: \n");	
		for (i = 0; i < m->size(); i++) {
			sum += localCSM[i];
			fprintf(out,"%s %7lf\n", m->symbol(i), localCSM[i]);
		}
		fprintf(out,"\nsum: %7lf\n", sum);
	}


}

/*
* prints PDB ATOM tags
*/
void printPDBATOM(Molecule* m,FILE* f,char** sym,double** pos){
	int i;
	for(i=0; i<m->size(); i++){
		fprintf(f,"ATOM  %5d %2s                %8.3lf%8.3lf%8.3lf                      %2s\n",
			i+1,sym[i],pos[i][0],pos[i][1],pos[i][2],sym[i]);
	}
}

/*
* prints PDB CONECT tags
*/
void printPDBCONNECT(Molecule* m,FILE* f){
	int i,j;
	for(i=0; i<m->size(); i++){
		fprintf(f,"CONECT%5d",i +1);
		for ( j=0;  j< m->valency(i) ; j++ ){
			if ((j>0) && (!(j%4)))
				fprintf(f,"\nCONECT%5d",i +1);
			fprintf(f,"%5d",m->adjacent(i,j) +1);
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

	printPDBATOM(m,out,m->symbols(),m->pos());

	printPDBCONNECT(m,out);

	fprintf(out,"ENDMDL\n");
	fprintf(out,"MODEL        2\n");

	printPDBATOM(m,out,m->symbols(),outAtoms);

	printPDBCONNECT(m,out);

	fprintf(out,"ENDMDL\n");

	// print results to screen


	if (writeOpenu)
		printf("SV* %.4lf *SV\n",fabs(csm));
	else
		printf( "%s: %.4lf\n",opName,fabs(csm));

	if (printNorm) {
		printf( "NORMALIZATION FACTOR: %7lf\n", m->norm());
		printf( "SCALING FACTOR OF SYMMETRIC STRUCTURE: %7lf\n", dMin);
		printf( "DIRECTIONAL COSINES: %lf %lf %lf\n", dir[0], dir[1], dir[2]);
		printf( "NUMBER OF EQUIVALENCE GROUPS: %d\n", m->groupNum());
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

	updateCoordinates(mol, m->pos());

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
		printf( "NORMALIZATION FACTOR: %7lf\n", m->norm());
		printf( "SCALING FACTOR OF SYMMETRIC STRUCTURE: %7lf\n", dMin);
		printf( "DIRECTIONAL COSINES: %lf %lf %lf\n", dir[0], dir[1], dir[2]);
		printf( "NUMBER OF EQUIVALENCE GROUPS: %d\n", m->groupNum());
	}

	if (printLocal) {
		double sum = 0;
		int i;
		fprintf(out,"\nLocal CSM: \n");	
		for (i = 0; i < m->size(); i++) {
			sum += localCSM[i];
			fprintf(out,"%s %7lf\n", m->symbol(i), localCSM[i]);
		}
		fprintf(out,"\nsum: %7lf\n", sum);
	}}

