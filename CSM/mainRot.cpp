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
#include "options.h"
#include "groupPermuter.h"
#include <openbabel/mol.h>
#include "babelAdapter.h"
#include <vector>
#include "Molecule.h"

#include "logging.h"
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <iomanip>
#include <sstream>

#include "math_wrappers.h"
#include "math_utils.h"

using namespace std;

#include "dvector.h"
#include "dmatrix.h"

#include "PrintOuts.h"
#include "drand48.h"
#include "calculations.h"

#define CSMFORMAT "CSM"

csm_options options;

void readPerm(FILE* permfile, int* perm, int size);
void readDir(FILE* dirFile, double* dir);

int mainWithOptions(); // Runs the main code after the options have been set

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
	if (options.type != CH) {
		return numPermutations(m, options.opOrder, options.type);
	} else {
		// CS
		int i;
		double numPerms = numPermutations(m, 2, SN);
		for (i = 2; i <= options.sn_max; i += 2) {
			numPerms += numPermutations(m, i, SN);
		}	
		return numPerms;		
	}
}

const char *getExtension(const char *fname) {
	return strrchr(fname,'.') + 1;
}

// ************************************************************
//       main program
// ************************************************************

/*
* main funciton - check valid parameters, parse molecule and call chirality Operation
*/
int main(int argc, char *argv[])
{
	cout << "argc: " << argc << endl;
	options = csm_options(argc, argv);

	int rc = mainWithOptions();
	return rc;
}

int mainWithOptions()
{
	init_logging();
	LOG(info) << "CSM starting up";

	int i;
	double csm, dMin;
	double **outAtoms;                 // output atom coordinates
	double dir[3] = { 0.0, 0.0, 0.0 };   // directional cosines
	int *perm = NULL;
	double *localCSM = NULL;

	if (options.logFileName != "")
		set_file_logging(options.logFileName);

	if ((options.findPerm && options.useperm) || (options.findPerm && options.useDir) || (options.useDir && options.useperm)) {
		LOG(fatal) << "-findperm, -useperm and -usedir are mutually exclusive";
		exit(1);
	} 

	// try to read molecule from infile
	Molecule* m;
	OBMol mol; 
	OperationType chMinType = CS;
	int chMinOrder = 2;
	
	if (options.useFormat) {
		// If a specific format is used, read molecule using that format
		if (boost::iequals(options.format, CSMFORMAT)) // Case-insensitive comparison
		{
			m = Molecule::create(options.inFile, stdout, options.ignoreSym && !options.useperm);
			if (m==NULL) exit(1);
			if (options.useMass)
			{
				m->fillAtomicMasses();
			}
		} else {
			mol = readMolecule(options.inFileName.c_str(), options.format, options.babelBond);
			m = Molecule::createFromOBMol(mol, options.ignoreSym && !options.useperm, options.useMass);
		}
   	} else {
		options.format = getExtension(options.inFileName.c_str());

		// if the extension is CSM - use csm
		if (boost::iequals(options.format, CSMFORMAT)) {
			m = Molecule::create(options.inFile, stdout, options.ignoreSym && !options.useperm);
			if (m==NULL) exit(1);
			if (options.useMass)
			{
				m->fillAtomicMasses();
			}
		} else {
			
			mol = readMolecule(options.inFileName.c_str(), "", options.babelBond);
			m = Molecule::createFromOBMol(mol, options.ignoreSym && !options.useperm, options.useMass);
		}
   	}

	if (options.babelTest) // Mol is ok - return 0
		return 0;

	if (!m){
		if (options.writeOpenu) {
			printf("ERR* Failed to read molecule from data file *ERR\n");
		}
		LOG(fatal) << "Failed to read molecule from data file";
		exit(1);
	}

	// strip unwanted atoms if needbe
	if ((options.ignoreHy || options.removeHy) && !options.useperm){
		char* removeList[] = {"H"," H"};
		Molecule* n = NULL;
		if (options.ignoreHy)
			n = m->stripAtoms(removeList,2,false);
		else //removeHy 
			n = m->stripAtoms(removeList,2,true);		
	
		mol.DeleteHydrogens();
		
	
		if (!n){
			if (options.writeOpenu) {
				printf("ERR* Failed while trying to strip unwanted atoms *ERR\n");
			}
			LOG(fatal) << "Failed while trying to strip unwanted atoms";
			exit(1);
		}
		delete m;
		m = n;
		// continue as per usual
	}

	if (!options.findPerm) {
		if (!options.useperm && !options.useDir) {
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
		if (options.timeOnly) { return 0; };
	}

	// allocate memory for outAtoms
	outAtoms = (double **)malloc(m->size() * sizeof(double*));
	for (i=0;i<m->size();i++)
		outAtoms[i] = (double *)malloc(3 * sizeof(double));
       
	perm = (int *)malloc(sizeof(int) * m->size());

	//normalize Molecule
	if (!m->normalizeMolecule(options.keepCenter)){
		if (options.writeOpenu) {
			printf("ERR* Failed to normalize atom positions: dimension of set of points = zero *ERR\n");
		}
		LOG(fatal) << "Failed to normalize atom positions: dimension of set of points = zero";
		exit(1);
	}
 
	if (options.useDir) readDir(options.dirfile, dir);

	if (options.useperm) {
		if (options.type == CH) {
			LOG(fatal) << "Chirality can't be given a permutation, run the specific csm operation instead";
			exit(1);
		}	
		readPerm(options.permfile, perm, m->size());
		runSinglePerm(m, outAtoms, perm, &csm, dir, &dMin, options.type);
	} else {
		if (options.type != CH) {
			// perform operation
			if (options.useDir) {
				findBestPermUsingDir(m, outAtoms, perm, &csm, dir, &dMin, options.type);
			}
			else if (options.findPerm) {
				findBestPerm(m, outAtoms, perm, &csm, dir, &dMin, options.type);
			} else {
				csmOperation(m, outAtoms, perm, &csm, dir, &dMin, options.type);
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
			options.opOrder = 2;
			if (options.useDir) {
				findBestPermUsingDir(m, outAtoms, perm, &csm, dir, &dMin, CS);				
			}
			else if (options.findPerm) {
				findBestPerm(m, outAtoms, perm, &csm, dir, &dMin, CS);				
			} else { 
				csmOperation(m, outAtoms, perm, &csm, dir, &dMin, CS);			
			}

			if (csm > MINDOUBLE) {							
				for (i = 2; i <= options.sn_max; i += 2) {
					options.opOrder = i;
					if (options.useDir) {
						findBestPermUsingDir(m, chOutAtoms, chPerm, &chCsm, chDir, &chdMin, SN);
					}
					else if (options.findPerm) {
						findBestPerm(m, chOutAtoms, chPerm, &chCsm, chDir, &chdMin, SN);
					} else {
						csmOperation(m, chOutAtoms, chPerm, &chCsm, chDir, &chdMin, SN);
					}				
					if (chCsm < csm) {
						int j;
						chMinType = SN;
						chMinOrder = options.opOrder;
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

	if (options.printLocal) {
		localCSM = (double *)malloc(sizeof(double) * m->size());
		if (options.type == CH) 
			options.opOrder = chMinOrder;
		computeLocalCSM(m, localCSM, perm, dir, options.type != CH ? options.type : chMinType);
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

	if (options.useFormat) {
		// If a specific format is used, read molecule using that format
		if (boost::iequals(options.format, CSMFORMAT)) // Case insensitive comparison 
		{
			printOutput(m, outAtoms, csm, dir, dMin, options.outFile, localCSM);
		} else {
			printOutputFormat(m, mol, outAtoms, csm, dir, dMin, options.outFile, options.outFileName.c_str(), localCSM);
		}
	} else {
		// if the extension is CSM - use csm
		if (strcasecmp(getExtension(options.inFileName.c_str()), CSMFORMAT) == 0) {
			printOutput(m, outAtoms, csm, dir, dMin, options.outFile, localCSM);
		} else {
			printOutputFormat(m, mol, outAtoms, csm, dir, dMin, options.outFile, options.outFileName.c_str(), localCSM);
		}
	}

	if (options.type == CH) {
		if (chMinType == CS) { 		
			fprintf(options.outFile, "\n MINIMUM CHIRALITY WAS FOUND IN CS\n\n");
		} else { 
			fprintf(options.outFile, "\n MINIMUM CHIRALITY WAS FOUND IN S%d\n\n", chMinOrder);
		}
	}	

	fprintf(options.outFile, "\n PERMUTATION:\n\n");
	for (i = 0; i < m->size(); i++) {
		fprintf(options.outFile, "%d ", perm[i] + 1);
	}
	fprintf(options.outFile, "\n");

	// housekeeping
	for (i=0;i<m->size();i++){
		free(outAtoms[i]);
	}
	free(outAtoms);
	delete m;
	free(perm);

	if (options.printLocal) free(localCSM);

	fclose(options.inFile);
	fclose(options.outFile);

	if (options.permfile != NULL)
		fclose(options.permfile);

	if (options.dirfile != NULL)
		fclose(options.dirfile);

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
			if (options.writeOpenu) {
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

