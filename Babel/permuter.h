/*
 * Author: shadi lahham
 *
 * generates permutations in the range 1 .. size
 * the permutations generated are limited to those containing
 * closed permutation groups of size groupSize
 *
 * permuter can generate permutations from a given integer array
 *
 */

#ifndef PERMUTER_H
#define PERMUTER_H

#define TRUE 1
#define FALSE 0

typedef struct cycle {

	// This cycle's size
	int _cycleSize;

	// The index into the divisor array
	int _divisor;

	// The actual cycle entries (_cycle[i]->_cycle[(i+1)%_cycleSize])
	int *_cycle;
} cycle;

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

permuter* createPermuter(int size,int groupSize, int addGroupsOfTwo);

int nextPermutation(permuter *p);

void resetPermuter(permuter *p);

void freePermuter(permuter *p);

#endif
