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
#include "drand48.h"
#include "calculations.h"

#define CSMFORMAT "CSM"

csm_options options;
csm_output results;

int mainWithOptions(); // Runs the main code after the options have been set
static void fill_output(Molecule *m, double **outAtoms, double csm, double *dir, double dMin, double *localCSM, 
	int chMinOrder, OperationType chMinType, int *perm);

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

	// try to read molecule from infile
	Molecule* m = options.molecule;
	OperationType chMinType = CS;
	int chMinOrder = 2;

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
 
	if (options.useDir)
	{
		dir[0] = options.dir[0];
		dir[1] = options.dir[1];
		dir[2] = options.dir[2];
	}

	if (options.useperm) {
		for (int i = 0; i<options.perm.size(); i++){
			perm[i] = options.perm[i];
		}
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
			// parameters dictionary change after the calculations
			// TODO: check if it's used in python after calculations to python
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

	fill_output(m, outAtoms, csm, dir, dMin, localCSM, chMinOrder, chMinType, perm);

	// housekeeping
	for (i=0;i<m->size();i++){
		free(outAtoms[i]);
	}
	free(outAtoms);
	delete m;
	free(perm);

	if (options.printLocal) free(localCSM);

	return 0;
}

void fill_output(Molecule *m, double **outAtoms, double csm, double *dir, double dMin, double *localCSM, int chMinOrder, 
	OperationType chMinType, int *perm)
{
	// All arrays that depend on the molecule's size
	results.molecule.atoms.clear();
	results.outAtoms.clear();
	results.localCSM.clear();
	results.perm.clear();
	for (int i = 0; i < m->size(); i++)
	{
		// Molecule
		python_atom atom;
		atom.symbol = m->symbol(i);
		atom.mass = m->mass(i);
		for (int j = 0; j < 3; j++)
			atom.pos.push_back(m->pos()[i][j]);
		for (int j = 0; j < m->valency(i); j++)
			atom.adjacent.push_back(m->adjacent(i, j));
		results.molecule.atoms.push_back(atom);

		// outAtoms
		std::vector<double> outAtom;
		for (int j = 0; j < 3; j++)
			outAtom.push_back(outAtoms[i][j]);
		results.outAtoms.push_back(outAtom);

		// localCSM
		if (localCSM)
			results.localCSM.push_back(localCSM[i]);

		// perm
		results.perm.push_back(perm[i]);
	}

	// Molecule equivalencyClasses
	results.molecule.equivalenceClasses.clear();
	int *group = new int[m->size()]; // No group is larger than the molecule - this is enough
	for (int i = 1; i <= m->groupNum(); i++)
	{
		int groupSize = m->getGroup(i, group);
		results.molecule.equivalenceClasses.push_back(vector<int>(group, group + groupSize));
	}
	delete[] group;


	// Other values
	results.norm = m->norm();
	results.numGroups = m->groupNum();
	results.csm = csm;

	results.dir.clear();
	for (int i = 0; i < 3; i++)
		results.dir.push_back(dir[i]);

	results.dMin = dMin;
	results.chMinOrder = chMinOrder;
	switch (chMinType) {
		case CN: results.chMinType = "CN"; break;
		case SN: results.chMinType = "SN"; break;
		case CS: results.chMinType = "CS"; break;
		case CI: results.chMinType = "CI"; break;
		case CH: results.chMinType = "CH"; break;
	}
}
