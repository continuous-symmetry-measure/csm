/*
 * Author: Amir Zait
 *
 * generates permutations in the range 1 .. size
 * the permutations generated are limited to those containing
 * cycles with length that is a divisor of the operation order
 *
 */

#include <math.h>
#include <assert.h>
#include "permuter.h"

using namespace std;

 /************************************************************************/
 /* Utility Functions                                                    */
 /************************************************************************/


 /**
 * Locate the next free index starting at a given index, or -1 if not found
 * TODO: consider different data structure
 * 
 */
 int Permuter::findFreeIndex(int startingAt) {
	 int i;	
	 for (i = startingAt; i < _size; i++) {
		 if (!_used[i])
			 return i;
	 }
	 return -1;
}


 /**
 * Find the divisors of the operation order and put them into the array
 */
 void Permuter::findDivisors() {
	 int i;
	 int maxDivisor = _operationOrder;
	 _numDivisors = 0;
	 _divisors = vector<int>(maxDivisor);
	 for (i = 2; i <= maxDivisor; i++) {
		 if ((_operationOrder % i) == 0) {
			 _divisors[_numDivisors++] = i;
		 }
	 }
}

/** 
* Cancel a cycle - return all indices participating in it to 
* the identity.
*/
void Permuter::cancelCycle(Cycle *c, bool clearUsed) {
	int i, ind;
	for (i = 0; i < c->_cycleSize; i++) {
		ind = c->_cycle[i];
		_index[ind] = ind;
		if (clearUsed)
			_used[ind] = false;
	}
}

/** 
* This function simply applies the given cycle - 
* Causes the permutation to contain it
*/
void Permuter::applyCycle(Cycle *c) {
	int i;
	for (i = 0; i < c->_cycleSize; i++) {		
		_index[c->_cycle[i]] = c->_cycle[(i + 1) % c->_cycleSize];
		_used[c->_cycle[i]] = true;
	}
}

int Permuter::initCycle(Cycle* c, int divisorIndex, int startingIndex, int ignoreUsed) {
	int i;	
	c->_cycleSize = _divisors[divisorIndex];
	c->_cycle = vector<int>(c->_cycleSize);
	c->_divisor = divisorIndex; 
	c->_cycle[0] = findFreeIndex(startingIndex);
	for (i = 1; i < c->_cycleSize; i++) {
		if (ignoreUsed) {
			c->_cycle[i] = startingIndex + i;
		} else {
			c->_cycle[i] = findFreeIndex(c->_cycle[i-1] + 1);
			if (c->_cycle[i] == -1) 
				return false;
		}
	}
	return true;
}

void Permuter::freeCycle(Cycle *c) {
	delete c;
}

/**
 * Clear the used array
 */
void Permuter::clearPerm() {
	int i;
	for (i = 0; i < _size; i++) {
		_used[i] = false;
		_index[i] = i;
	}
}

/** 
 * Try to advance a cycle, return false upon failure
 */
int Permuter::advanceCycle(int cycleIndex) {
	Cycle* c = _cycles[cycleIndex];
	int k = c->_cycleSize - 1; 
	cancelCycle(c, false);
	while (k >= 0) {
		_used[c->_cycle[k]] = false;
		// If we can still advance this index - do it. 
		if (c->_cycle[k] != _size - 1) {
			int j = k;
			int freeIndex;
			while (j < c->_cycleSize) {
				// Make sure we are actually advancing, and remember that the first index 
				// in the cycle is always the smallest
				if (k == j) {
					freeIndex = findFreeIndex(c->_cycle[k] + 1); 
				} else {
					freeIndex = findFreeIndex(c->_cycle[0] + 1); 
				} 
				if (freeIndex == -1) 
					break;
				c->_cycle[j] = freeIndex;
				_used[freeIndex] = true;
				j++;
			}
			// if we advanced all cycles
			if (j == c->_cycleSize) {
				applyCycle(c);
				return true;
			} else {
				// unroll
				int i;
				for (i = k; i < j; i++) {
					_used[c->_cycle[i]] = false;
				}
			}
		}		
		k--;
	}
	return false;
}

/************************************************************************/
/* Interface Functions                                                                     */
/************************************************************************/

Permuter::Permuter(int size, int groupSize, int addGroupsOfTwo) : _used(size), _index(size)
{	
	_size = size;
	_operationOrder = groupSize;
	// The cycles can only be of length 1 or operation order. 
	//findDivisors(p);
	if (addGroupsOfTwo && groupSize > 2) {
		_divisors = { 2, _operationOrder };
		//_divisors[0] = 2;
		//_divisors[1] = _operationOrder;
		_numDivisors = 2;
	} else {
		_divisors = { _operationOrder };
		//_divisors[0] = _operationOrder;
		_numDivisors = 1;
	}

	_maxCycles = _size / _divisors[0];
	_cycles = vector<Cycle *>(_maxCycles);	
	_numCycles = 0;
	reset();
}

bool Permuter::next() {
	if (_firstPermutation) {
		_firstPermutation = false;
		return true;
	} else {
		int k = _numCycles - 1;
		while (k >= 0) {
			// try to advance current cycle
			if (advanceCycle(k)) {
				int curSize = _cycles[k]->_cycleSize;
				int curStart = _cycles[k]->_cycle[0] + 1;
				int i;
				// if successful - try to initialize cycles following it.
				for (i = k + 1; i < _numCycles; i++) {
					int sameSize = (_cycles[i]->_cycleSize == curSize);
					// This cycle has already been canceled, no need to re-cancel
					if (initCycle(_cycles[i],_cycles[i]->_divisor, sameSize ? curStart : 0, false)){
						applyCycle(_cycles[i]);
					} else {
						int j;
						for (j = k + 1; j < i; j++) {
							cancelCycle(_cycles[j], true);
						}
						break;
					}
					curSize = _cycles[i]->_cycleSize;
					curStart = _cycles[i]->_cycle[0] + 1;
				}
				if (i == _numCycles) {
					return true;
				} else {
					cancelCycle(_cycles[k], true);
				}
			} 
			k--;
		}
		// If we didn't succeed in advancing any cycle, begin from last one, again, and
		// try to increase its size
		k = _numCycles - 1;
		while (k >= 0) {
			if (_cycles[k]->_divisor + 1 < _numDivisors) {
				int i;
				int numUsed = 0;
				int newDivisorIndex = _cycles[k]->_divisor + 1;
				for (i = 0; i < k; i++) {		
					numUsed += _divisors[_cycles[i]->_divisor];
				}
				for (i = k; i < _numCycles; i++) {
					numUsed += _divisors[newDivisorIndex];
				}
				// if possible - make this circle and its followers larger by one				
				if (numUsed <= _size) {
					int numUsed = 0;
					int newDivisorIndex = _cycles[k]->_divisor + 1;
					int i;
					clearPerm();
					// up to the cycle who's size is changing - just initialize permutation
					for (i = 0; i < k; i++) {		
						initCycle(_cycles[i], _cycles[i]->_divisor, numUsed, true);
						numUsed += _divisors[_cycles[i]->_divisor];
						applyCycle(_cycles[i]);
					}
					// from this cycle and forward - all cycles are the same (new) size;
					for (i = k; i < _numCycles; i++) {
						initCycle(_cycles[i], newDivisorIndex, numUsed, true);
						numUsed += _divisors[newDivisorIndex];
						applyCycle(_cycles[i]);
					}
					_numUsed = numUsed;
					return true;
				}
			}
			k--;
		}


		// we've reached the end 
		if (_numCycles >= _maxCycles) {
			// if no more cycles can be added
			return false;
		} else {
			// add a cycle, resize all cycles to smallest size			
			int i;
			clearPerm();
			_cycles[_numCycles] = new Cycle();
			// Don't initialize this, initCycle below will take care of it
			//_cycles[_numCycles]->_cycle = (int *)malloc(_operationOrder * sizeof(int));
			_numCycles++;
			for (i = 0; i < _numCycles; i++) {
				initCycle(_cycles[i], 0, i * _divisors[0], true);
				applyCycle(_cycles[i]);
			}
			_numUsed = _numCycles * _divisors[0];
		}
		return true;
	}
}

void Permuter::reset() 
{
	int i;
	_firstPermutation = true;
	for (i = 0; i < _numCycles; i++) {
		freeCycle(_cycles[i]);
		_cycles[i] = nullptr;
	}
	_numCycles = 0;
	for (i = 0; i < _size; i++) {
		_used[i] = 0;
		_index[i] = i;
	}
	_numUsed = 0;
}

Permuter::~Permuter()
{
	reset();
}

const int Permuter::operator[](int index) const
{
	assert(index >= 0 && index < _size);
	return _index[index];
}
