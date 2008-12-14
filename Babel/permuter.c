/*
 * Author: Amir Zait
 *
 * generates permutations in the range 1 .. size
 * the permutations generated are limited to those containing
 * cycles with length that is a divisor of the operation order
 *
 */

#include <stdlib.h>  // for NULL
#include <stdio.h>   // for printf
#include <math.h>
#include "permuter.h"

 /************************************************************************/
 /* Utility Functions                                                    */
 /************************************************************************/


 /**
 * Locate the next free index starting at a given index, or -1 if not found
 * TODO: consider different data structure
 * 
 */
 int findFreeIndex(permuter *p, int startingAt) {
	 int i;	
	 for (i = startingAt; i < p->_size; i++) {
		 if (!p->_used[i])
			 return i;
	 }
	 return -1;
}


 /**
 * Find the divisors of the operation order and put them into the array
 */
 void findDivisors(permuter *p) {
	 int i;
	 int maxDivisor = p->_operationOrder;
	 p->_numDivisors = 0;
	 p->_divisors = (int *)malloc(maxDivisor * sizeof(int));
	 for (i = 2; i <= maxDivisor; i++) {
		 if ((p->_operationOrder % i) == 0) {
			 p->_divisors[p->_numDivisors++] = i;
		 }
	 }
}

/** 
* Cancel a cycle - return all indices participating in it to 
* the identity.
*/
void cancelCycle(permuter *p, cycle *c, int clearUsed) {
	int i, ind;
	for (i = 0; i < c->_cycleSize; i++) {
		ind = c->_cycle[i];
		p->_index[ind] = ind;
		if (clearUsed)
			p->_used[ind] = FALSE;
	}
}

/** 
* This function simply applies the given cycle - 
* Causes the permutation to contain it
*/
void applyCycle(permuter *p, cycle *c) {
	int i;
	for (i = 0; i < c->_cycleSize; i++) {		
		p->_index[c->_cycle[i]] = c->_cycle[(i + 1) % c->_cycleSize];
		p->_used[c->_cycle[i]] = TRUE;
	}
}

int initCycle(permuter *p, cycle* c, int divisorIndex, int startingIndex, int ignoreUsed) {
	int i;	
	c->_cycleSize = p->_divisors[divisorIndex];
	if (c->_cycle != NULL && c->_divisor != divisorIndex) {
		free(c->_cycle);
		c->_cycle = (int *)malloc(sizeof(int) * c->_cycleSize);
	} else if (c->_cycle == NULL) {
		c->_cycle = (int *)malloc(sizeof(int) * c->_cycleSize);
	}
	c->_divisor = divisorIndex; 
	c->_cycle[0] = findFreeIndex(p, startingIndex);
	for (i = 1; i < c->_cycleSize; i++) {
		if (ignoreUsed) {
			c->_cycle[i] = startingIndex + i;
		} else {
			c->_cycle[i] = findFreeIndex(p, c->_cycle[i-1] + 1);
			if (c->_cycle[i] == -1) 
				return FALSE;
		}
	}
	return TRUE;
}

void freeCycle(cycle *c) {
	free(c->_cycle);
	free(c);
}

/**
 * Clear the used array
 */
void clearPerm(permuter *p) {
	int i;
	for (i = 0; i < p->_size; i++) {
		p->_used[i] = FALSE;
		p->_index[i] = i;
	}
}

/** 
 * Try to advance a cycle, return FALSE upon failure
 */
int advanceCycle(permuter *p, int cycleIndex) {
	cycle* c = p->_cycles[cycleIndex];
	int k = c->_cycleSize - 1; 
	cancelCycle(p, c, FALSE);
	while (k >= 0) {
		p->_used[c->_cycle[k]] = FALSE;
		// If we can still advance this index - do it. 
		if (c->_cycle[k] != p->_size - 1) {
			int j = k;
			int freeIndex;
			while (j < c->_cycleSize) {
				// Make sure we are actually advancing, and remember that the first index 
				// in the cycle is always the smallest
				if (k == j) {
					freeIndex = findFreeIndex(p, c->_cycle[k] + 1); 
				} else {
					freeIndex = findFreeIndex(p, c->_cycle[0] + 1); 
				} 
				if (freeIndex == -1) 
					break;
				c->_cycle[j] = freeIndex;
				p->_used[freeIndex] = TRUE;
				j++;
			}
			// if we advanced all cycles
			if (j == c->_cycleSize) {
				applyCycle(p, c);
				return TRUE;
			} else {
				// unroll
				int i;
				for (i = k; i < j; i++) {
					p->_used[c->_cycle[i]] = FALSE;
				}
			}
		}		
		k--;
	}
	return FALSE;
}

/************************************************************************/
/* Interface Functions                                                                     */
/************************************************************************/

permuter* createPermuter(int size,int groupSize, int addGroupsOfTwo) {	
	permuter* p = (permuter*)malloc(sizeof(permuter));
	p->_size = size;
	p->_operationOrder = groupSize;
	// The cycles can only be of length 1 or operation order. 
	//findDivisors(p);
	if (addGroupsOfTwo && groupSize > 2) {
		p->_divisors = (int *)malloc(sizeof(int) * 2);
		p->_divisors[0] = 2;
		p->_divisors[1] = p->_operationOrder;
		p->_numDivisors = 2;
	} else {
		p->_divisors = (int *)malloc(sizeof(int) * 1);
		p->_divisors[0] = p->_operationOrder;
		p->_numDivisors = 1;
	}

	p->_maxCycles = p->_size / p->_divisors[0];
	p->_cycles = (cycle**)malloc(sizeof(cycle*) * p->_maxCycles);	
	p->_index = (int *)malloc(sizeof(int) * size);
	p->_used = (int *)malloc(sizeof(int) * size);
	p->_numCycles = 0;
	resetPermuter(p);
	return p;
}

int nextPermutation(permuter *p) {
	if (p->_firstPermutation) {
		p->_firstPermutation = FALSE;
		return TRUE;
	} else {
		int k = p->_numCycles - 1;
		while (k >= 0) {
			// try to advance current cycle
			if (advanceCycle(p, k)) {
				int curSize = p->_cycles[k]->_cycleSize;
				int curStart = p->_cycles[k]->_cycle[0] + 1;
				int i;
				// if successful - try to initialize cycles following it.
				for (i = k + 1; i < p->_numCycles; i++) {
					int sameSize = (p->_cycles[i]->_cycleSize == curSize);
					// This cycle has already been canceled, no need to re-cancel
					if (initCycle(p, p->_cycles[i],p->_cycles[i]->_divisor, sameSize ? curStart : 0, FALSE)){
						applyCycle(p, p->_cycles[i]);
					} else {
						int j;
						for (j = k + 1; j < i; j++) {
							cancelCycle(p, p->_cycles[j], TRUE);
						}
						break;
					}
					curSize = p->_cycles[i]->_cycleSize;
					curStart = p->_cycles[i]->_cycle[0] + 1;
				}
				if (i == p->_numCycles) {
					return TRUE;
				} else {
					cancelCycle(p, p->_cycles[k], TRUE);
				}
			} 
			k--;
		}
		// If we didn't succeed in advancing any cycle, begin from last one, again, and
		// try to increase its size
		k = p->_numCycles - 1;
		while (k >= 0) {
			if (p->_cycles[k]->_divisor + 1 < p->_numDivisors) {
				int i;
				int numUsed = 0;
				int newDivisorIndex = p->_cycles[k]->_divisor + 1;
				for (i = 0; i < k; i++) {		
					numUsed += p->_divisors[p->_cycles[i]->_divisor];
				}
				for (i = k; i < p->_numCycles; i++) {
					numUsed += p->_divisors[newDivisorIndex];
				}
				// if possible - make this circle and its followers larger by one				
				if (numUsed <= p->_size) {
					int numUsed = 0;
					int newDivisorIndex = p->_cycles[k]->_divisor + 1;
					int i;
					clearPerm(p);
					// up to the cycle who's size is changing - just initialize permutation
					for (i = 0; i < k; i++) {		
						initCycle(p, p->_cycles[i], p->_cycles[i]->_divisor, numUsed, TRUE);
						numUsed += p->_divisors[p->_cycles[i]->_divisor];
						applyCycle(p, p->_cycles[i]);
					}
					// from this cycle and forward - all cycles are the same (new) size;
					for (i = k; i < p->_numCycles; i++) {
						initCycle(p, p->_cycles[i], newDivisorIndex, numUsed, TRUE);
						numUsed += p->_divisors[newDivisorIndex];
						applyCycle(p, p->_cycles[i]);
					}
					p->_numUsed = numUsed;
					return TRUE;
				}
			}
			k--;
		}


		// we've reached the end 
		if (p->_numCycles >= p->_maxCycles) {
			// if no more cycles can be added
			return FALSE;
		} else {
			// add a cycle, resize all cycles to smallest size			
			int i;
			clearPerm(p);
			p->_cycles[p->_numCycles] = (cycle*)malloc(sizeof(cycle));
			p->_cycles[p->_numCycles]->_cycle = (int *)malloc(p->_operationOrder * sizeof(int));
			p->_numCycles++;
			for (i = 0; i < p->_numCycles; i++) {
				initCycle(p, p->_cycles[i], 0, i * p->_divisors[0], TRUE);
				applyCycle(p, p->_cycles[i]);
			}
			p->_numUsed = p->_numCycles * p->_divisors[0];
		}
		return TRUE;
	}
}

void resetPermuter(permuter *p) {
	int i;
	p->_firstPermutation = TRUE;
	for (i = 0; i < p->_numCycles; i++) {
		freeCycle(p->_cycles[i]);
		p->_cycles[i] = NULL;
	}
	p->_numCycles = 0;
	for (i = 0; i < p->_size; i++) {
		p->_used[i] = 0;
		p->_index[i] = i;
	}
	p->_numUsed = 0;
}

void freePermuter(permuter *p) {
	resetPermuter(p);
	free(p->_cycles);
	free(p->_index);
	free(p->_used);
	free(p->_divisors);
	free(p);
}

