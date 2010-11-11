#include <stdio.h>
#include <python.h>

using namespace std;

extern "C" {
	#include "Molecule.h"
}
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
extern FILE* permfile = NULL;

extern char opName[100];

/**
 * The arguments are in the following form:
 * args = (positions, elements, csm_type, args)
 */
PyObject *computeCsm(PyObject* self, PyObject* args) {
	PyObject *coords = NULL;
	PyObject *elems = NULL;	
	char *csmType = NULL;
	PyObject *optionList = NULL;
	Molecule *mol = NULL;

	// Read the arguments
	if ( ! PyArg_ParseTuple(args, "(OOsO)", &coords, &elems, &csmType, &optionList) ) {
		printf("Could not unparse objects\n");
        return NULL;
	}	

	// Construct molecule
	mol = PyMol2Mol(coords, elems);

	// Return the result as a tuple
	// (csm, localCSM, );
	return PyMain(mol, csmType, optionList);
}

/** 
 * The molecular data will be given as two lists:
 * tuples of (x,y,z)
 * list of element names inside tuples (name) 
 */
Molecule* PyMol2Mol(PyObject *coords, PyObject *elements) {
  int numAtoms = PyList_Size(coords);  
  Molecule *mol = allocateMolecule(numAtoms);
    
  for (int i = 0; i < numAtoms; i++) { 
    PyObject* coords = PyList_GetItem(coords, i);
    PyObject* elem = PyList_GetItem(element, i);
    PyArg_ParseTuple("fff", &(mol->_pos[i][0]), &(mol->_pos[i][1]), &(mol->_pos[i][2]));
    PyArg_ParseTuple("s", ,&(m->_symbol[i]));        
	
	// So far - we don't know how to extract valency and bond structure...
	mol->_valency[i] = 0;
	mol->_adjacent[i] = (int*)malloc(0);
  }
  
  initSimilarity(mol,DEPTH_ITERATIONS);
  return mol;
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
		opOrder = atoi(argv[1] + 1);
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

	for ( i=0;  i< numItems ;  i++ ){
		PyObject* elem = PyList_GetItem(element, i);
		char *arg = NULL;
		PyArg_ParseTuple("s", ,&(arg));  
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

		else if (strcmp(arg,"-ignoreSym" ) == 0 )
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
			usage(argv[0]);
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
		}
	}
	if (writeOpenu) {
		useFormat = TRUE;
		format = strdup("PDB");		
	}

	return TRUE;
}


/*
* main funciton - check valid parameters, parse molecule and call chirality Operation
*/
PyObject *PyMain(Molecule *mol, char *csmType, PyObject *optionList){

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
	
		mol.DeleteHydrogens();		
	
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

	if (printLocal) {	
		localCSM = (double *)malloc(sizeof(double) * m->_size);
		computeLocalCSM(m,localCSM, perm, dir,  type != CH ? type : chMinType);
	}

	
	// De-normalize
	for (i = 0; i < m->_size; i++) { 
		m->_pos[i][0] *= m->_norm;
		m->_pos[i][1] *= m->_norm;
		m->_pos[i][2] *= m->_norm;
		outAtoms[i][0] *= m->_norm;
		outAtoms[i][1] *= m->_norm;
		outAtoms[i][2] *= m->_norm;
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
	return Py_BuildValue("(d)", csm);
}