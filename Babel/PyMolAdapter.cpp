#include <stdio.h>
#include <python.h>

extern "C" {
	#include "Molecule.h"
}
typedef enum {
	CN,
	SN, 
	CS,
	CI, 
	CH
} OperationType;

#define MAXDOUBLE  100000000.0
#define MINDOUBLE  1e-8
#define TRUE 1
#define FALSE 0
#define APPROX_RUN_PER_SEC 8e4

extern Molecule * allocateMolecule(int size);
extern void initSimilarity(Molecule *m,int depth);		
extern void replaceSymbols(Molecule* m);

// function declarations
extern void csmOperation(Molecule* m, double** outAtoms, int *optimalPerm, double* csm, double* dir, double* dMin, OperationType type);
extern void runSinglePerm(Molecule* m, double** outAtoms, int* perm, double* csm, double* dir, double* dMin, OperationType type);
extern void findBestPerm(Molecule* m, double** outAtoms, int* optimalPerm, double* csm, double* dir, double* dMin, OperationType type);
extern void findSymmetryDirection(Molecule *m, double  ***dirs, int *n_dirs, OperationType type);
extern void estimatePerm(Molecule* m, int *perm, double *dir, OperationType type);
extern void readPerm(FILE* permfile, int* perm, int size);
extern double createSymmetricStructure(Molecule* m, double **outAtom, int *perm, double *dir, OperationType type, double dMin);
extern double computeLocalCSM(Molecule* m, double *localCSM, int *perm, double *dir, OperationType type);
extern double totalNumPermutations(Molecule *m);
extern void normalize(double **coords, Molecule *m);

// global options
extern int ignoreHy;
extern int removeHy;
extern int ignoreSym;
extern int writeOpenu;
extern OperationType type;
extern OperationType chMinType;
extern int opOrder;
extern int useperm;
extern int findPerm; 
extern int useMass;
extern int limitRun;
extern char *format;
extern int babelBond;
extern int timeOnly;
extern int sn_max;
extern int anneal;
extern int detectOutliers;
extern double A;
extern int babelTest;
extern int printNorm;
extern int printLocal;

// file pointers
extern FILE* permfile;

extern char opName[100];

Molecule* PyMol2Mol(PyObject *coords, PyObject *elements);
PyObject *PyMain(Molecule *mol, char *csmType, PyObject *optionList);
Molecule* PyMol2Mol(PyObject *coords, PyObject *elements, PyObject *bonds);

/**
 * The arguments are in the following form:
 * args = (positions, elements, bonds, csm_type, args)
 */
PyObject *computeCsm(PyObject* self, PyObject* args) {
	PyObject *coords = NULL;
	PyObject *elems = NULL;
	PyObject *bonds = NULL;
	char *csmType = NULL;
	PyObject *optionList = NULL;
	Molecule *mol = NULL;

	// Read the arguments
	if ( ! PyArg_ParseTuple(args, "OOOsO", &coords, &elems, &bonds, &csmType, &optionList) ) {
		printf("Could not unparse objects\n");
        return NULL;
	}	

	// Construct molecule
	mol = PyMol2Mol(coords, elems, bonds);

	// Return the result as a tuple
	// (csm, localCSM, );
	return PyMain(mol, csmType, optionList);	
}

/** 
 * The molecular data will be given as two lists:
 * tuples of (x,y,z)
 * list of element names inside tuples (name) 
 */
Molecule* PyMol2Mol(PyObject *coords, PyObject *elements, PyObject *bonds) {
  int numAtoms = PyList_Size(coords);  
  Molecule *mol = allocateMolecule(numAtoms);
      
  for (int i = 0; i < numAtoms; i++) { 
    PyObject *coord = PyList_GetItem(coords, i);
    PyObject *elem = PyList_GetItem(elements, i);	
	PyObject *atomBonds = PyList_GetItem(bonds, i);

	mol->_pos[i][0] = PyFloat_AsDouble(PyTuple_GetItem(coord, 0));
	mol->_pos[i][1] = PyFloat_AsDouble(PyTuple_GetItem(coord, 1));
	mol->_pos[i][2] = PyFloat_AsDouble(PyTuple_GetItem(coord, 2));	
	mol->_symbol[i] = strdup(PyString_AsString(elem));
	
	// So far - we don't know how to extract valency and bond structure...
	mol->_valency[i] = PyList_Size(atomBonds);
	mol->_adjacent[i] = (int*)malloc(mol->_valency[i] * sizeof(int));
	for (int j = 0; j < mol->_valency[i]; j++) {
		mol->_adjacent[i][j] = PyInt_AsLong(PyList_GetItem(atomBonds, j));
	}
  }
    
  initSimilarity(mol,DEPTH_ITERATIONS);
  return mol;
}

void PyUsage(char *op) {
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
	printf("-anneal		   - Try to anneal the result\n");
	printf("-babelbond	   - Let openbabel compute bonding\n");
	printf("-useMass	   - Use the atomic masses to define center of mass\n");
	printf("-timeOnly	   - Only print the time and exit\n");
	printf("-babelTest	   - Test if the molecule is legal or not\n");
	printf("-sn_max	<max n> - The maximal sn to try, relevant only for chirality\n");
	printf("-printNorm		- Print the normalization factor as well\n");
	printf("-printlocal		- Print the local CSM (csm for each atom) in the output file\n");
	printf("-help - print this help file\n");
}



/*
* parses the command line parameters
*/
int parsePyInput(char *csm_type, PyObject *args) {

	if (!strcmp(csm_type,"cs")) {
		type = CS;
		opOrder = 2;
		sprintf(opName, "MIRROR SYMMETRY");			
	} else if (!strcmp(csm_type, "ci")) {
		type = CI;
		opOrder = 2;
		sprintf(opName, "INVERSION (S2)");	
	} else if (!strcmp(csm_type, "ch")) {
		type = CH;
		opOrder = 2;
		sprintf(opName, "CHIRALITY");	
	} else if (csm_type[0] == 'c') {
		type = CN;
		opOrder = atoi(csm_type + 1);
		sprintf(opName, "C%d SYMMETRY",opOrder);	
	} else if (csm_type[0] == 's') {
		type = SN;
		opOrder = atoi(csm_type + 1);
		if (opOrder % 2 != 0) {
			printf("ERROR - Only Even values of n are allowed\n");
			return FALSE;
		}
		sprintf(opName, "S%d SYMMETRY",opOrder);	
	}	
	
	// get commandline flags
	int numItems = PyList_Size(args);  
	int nextIsPermFile = FALSE;
	int nextIsMaxSn = FALSE;

	for (int i=0;  i< numItems ;  i++ ){
		PyObject* elem = PyList_GetItem(args, i);
		char *arg = PyString_AsString(elem);		
		if (nextIsPermFile) {
			char* permfileName = arg;
			if ((permfile = fopen(permfileName, "rt")) == NULL){
				if (writeOpenu) {
					printf("ERR* Failed to open perm file %s for reading *ERR\n", permfileName);
				} else {
					printf("Failed to open perm file %s for reading\n", permfileName);
				}
				return FALSE;
			}
			nextIsPermFile = FALSE;
		} else if (nextIsMaxSn) { 
			sn_max = atoi(arg);
			nextIsMaxSn = FALSE;
		} else if (strcmp(arg,"-sn_max" ) == 0) {
			if (type != CH) { 
				printf("This option only applies to chirality\n");
				return FALSE;
			}
			nextIsMaxSn = TRUE;	
		} else if (strcmp(arg,"-ignoreHy" ) == 0 )
			ignoreHy = TRUE;
		else if (strcmp(arg,"-removeHy" ) == 0 )
			removeHy = TRUE;

		else if (strcmp(arg,"-ignoreSym" ) == 0 ) {
			ignoreSym = TRUE;
		} else if (strcmp(arg, "-nolimit") == 0) { 
			limitRun = FALSE;  
		} else if (strcmp(arg, "-useperm") == 0) {
			useperm = TRUE;
			nextIsPermFile = TRUE;
		} else if (strcmp(arg, "-babelbond") == 0) {
			babelBond = TRUE;
		} else if (strcmp(arg, "-useMass") == 0) { 
			useMass = TRUE;					
		} else if (strcmp(arg, "-timeonly") == 0) {
			timeOnly = TRUE;
		} else if  (strcmp(arg, "-printNorm") == 0) {
			printNorm = TRUE;
		} else if (strcmp(arg, "-help") == 0) {
			PyUsage(csm_type);
			return FALSE;
		} else if (strcmp(arg, "-findperm") == 0) { 
			findPerm = TRUE;
		} else if (strcmp(arg, "-detectOutliers") == 0) {
			detectOutliers = TRUE;
		} else if (strcmp(arg, "-anneal") == 0) {
			anneal = TRUE;
		} else if (strcmp(arg, "-babelTest") == 0) { 
			babelTest = TRUE;
		} else if (strcmp(arg, "-printlocal") == 0) { 
			printLocal = TRUE;
		} else if (strcmp(arg, "-approx") == 0) { 
			findPerm = TRUE;
			detectOutliers = TRUE;
		}
	}
	if (writeOpenu) {		
		format = strdup("PDB");		
	}

	return TRUE;
}


/*
* main funciton - check valid parameters, parse molecule and call chirality Operation
*/
PyObject *PyMain(Molecule *m, char *csmType, PyObject *optionList) {

	int i;
	double csm, dMin;
	double **outAtoms;                 // output atom coordinates
	double dir[3] = {0.0, 0.0, 0.0};   // directional cosines
	int *perm = NULL;	
	double *localCSM = NULL;		 

	if (parsePyInput(csmType, optionList) == FALSE) { 
		printf("Could not process input options objects\n");
        return NULL;
	}

	if (findPerm && useperm) { 
		printf("-findperm and -useperm can't be used together...");
		return NULL;
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
			return NULL;
		}
		freeMolecule(m);
		m = n;
		// continue as per usual
	}

	if (!findPerm) {
		if (!useperm) {
			double time = 1.0*totalNumPermutations(m) / 3600 / APPROX_RUN_PER_SEC;
			if (time != time) {
				// time is NaN
				time = MAXDOUBLE;
			}
			printf("Going to enumerate over %5.2f permutations\n", totalNumPermutations(m));
			printf("Entire run should take approx. %5.2f hours on a 2.0Ghz Computer\n", time  );
		} else {
			printf("Going to enumerate over %5.2f permutations\n", 1.0);
			printf("Entire run should take approx. %5.2f hours on a 2.0Ghz Computer\n", 1.0 / 3600 / APPROX_RUN_PER_SEC );
		}
		if (timeOnly) { return NULL; };
	}	

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
		return NULL;
	}

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
			if (findPerm) { 
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

			chOutAtoms = (double **)malloc(m->_size * sizeof(double*));
			for (i=0;i<m->_size;i++)
				chOutAtoms[i] = (double *)malloc(3 * sizeof(double));
       
			chPerm = (int *)malloc(sizeof(int) * m->_size);
			chMinType = CS;						
			opOrder = 2;
			if (findPerm) { 
				findBestPerm(m, outAtoms, perm, &csm, dir, &dMin, CS);				
			} else { 
				csmOperation(m, outAtoms, perm, &csm, dir, &dMin, CS);			
			}

			if (csm > MINDOUBLE) {							
				for (i = 2; i <= sn_max; i+=2) {
					opOrder = i;
					if (findPerm) { 
						findBestPerm(m, chOutAtoms, chPerm, &chCsm, chDir, &chdMin, SN);
					} else {
						csmOperation(m, chOutAtoms, chPerm, &chCsm, chDir, &chdMin, SN);
					}
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

	printLocal = TRUE;
	if (printLocal) {	
		localCSM = (double *)malloc(sizeof(double) * m->_size);
		computeLocalCSM(m,localCSM, perm, dir,  type != CH ? type : chMinType);
	}

	normalize(outAtoms,m);

	
	// De-normalize
	for (i = 0; i < m->_size; i++) { 
		m->_pos[i][0] *= m->_norm;
		m->_pos[i][1] *= m->_norm;
		m->_pos[i][2] *= m->_norm;
		outAtoms[i][0] *= m->_norm;
		outAtoms[i][1] *= m->_norm;
		outAtoms[i][2] *= m->_norm;
	}	

	// Create a list of the new items
	PyObject *localCSMList = PyList_New(m->_size);
	for (i = 0; i < m->_size; i++) {
		PyList_SetItem(localCSMList, i, Py_BuildValue("f", localCSM[i]));
	}

	// create a list for the axis coordinates
	PyObject *pyDir = PyTuple_New(3);
	for (i = 0; i < 3; i++) {
		PyTuple_SetItem(pyDir, i, Py_BuildValue("f", dir[i]));
	}

	// Create a list for the outgoing coordinates
	PyObject *pyOutAtoms = PyList_New(m->_size);
	for (i = 0; i < m->_size; i++) {
		PyList_SetItem(pyOutAtoms, i, Py_BuildValue("(fff)", outAtoms[i][0],outAtoms[i][1], outAtoms[i][2]));
	}

	// housekeeping
	for (i=0;i<m->_size;i++){
		free(outAtoms[i]);
	}

	free(outAtoms);
	freeMolecule(m);
	free(perm);

	if (printLocal) free(localCSM);

	if (permfile != NULL)
		fclose(permfile);	

	// Construct the object to return
	return Py_BuildValue("(dOOO)", csm, localCSMList,pyDir,pyOutAtoms);
}

static PyMethodDef csm_methods[] = {
	{"computeCsm", computeCsm, METH_VARARGS, "computeCsm()"},
	{NULL, NULL}
};

PyMODINIT_FUNC
initcsm(void)
{
	Py_InitModule("csm", csm_methods);
}
