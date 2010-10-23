/*
 * Author: Shadi Lahham, modified by Amir Zait
 *
 * Enumerates over all permutations in the range 1 .. N, 
 * which are composed only of cycles of given sizes. 
 */

#ifndef PERMUTER_H
#define PERMUTER_H

#define TRUE 1
#define FALSE 0

// A structure representing a cycle in the permutation
typedef struct cycle {

	// This cycle's size
	int _cycleSize;

	// The index into the divisor array
	int _divisor;

	// The actual cycle entries (_cycle[i]->_cycle[(i+1)%_cycleSize])
	int *_cycle;
} cycle;

// The main structure for the permuter
typedef struct permuter {

	// The size of the permutation
	int _size;

	// The divisors
	int *_divisors;

	// The number of existing divisors 
	int _numDivisors;

	// The maximum number of cycles
	int _maxCycles;

	// The current number of cycles
	int _numCycles;
	
	// The indices used
	int *_used;

	// The number of indices participating in cycles
	int _numUsed;

	// The current cycle structure of the permutation
	cycle** _cycles;

	// Is this the first permutation
	int _firstPermutation;

	// The current permutation
	int *_index;

	// operation order
	int _operationOrder;
} permuter;

/** 
 * Initialize a new permuter
 * 
 * @param size The size of the permutations to generate (1 .. size)
 * @param groupSize The size of allowed cycles (in addition to cycles of size 1)
 * @param addGroupsOfTwo Should cycles of size 2 also be allowed
 * 
 * @return The new permuter object
 */
permuter* createPermuter(int size,int groupSize, int addGroupsOfTwo);

/**
 * Generate the next permutation
 * 
 * @param p The permuter object
 * 
 * @return TRUE if we have not yet reached the end of the enumerator, FALSE otherwise
 */
int nextPermutation(permuter *p);

/** 
 * Reset the permuter to re-initialize enumeration
 */
void resetPermuter(permuter *p);

/** 
 * Free the permuter
 */
void freePermuter(permuter *p);

#endif
