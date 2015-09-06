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