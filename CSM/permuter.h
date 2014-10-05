/*
 * Author: Shadi Lahham, modified by Amir Zait
 *
 * Enumerates over all permutations in the range 1 .. N, 
 * which are composed only of cycles of given sizes. 
 *
 * Converted to C++ by Itay Zandbank
 */

#ifndef PERMUTER_H
#define PERMUTER_H

#include <vector>

class Permuter
{
private:
	struct Cycle
	{
		// This cycle's size
		int _cycleSize;

		// The index into the divisor array
		int _divisor;

		// The actual cycle entries (_cycle[i]->_cycle[(i+1)%_cycleSize])
		std::vector<int> _cycle;

		/* C++ Conversion Note: Cycle wasn't converted into an actual C++ object because the code is a bit messy,
		* Cycle and Permuter logic is intertwined. Since this class is private and is just an implementation detail,
		* I decided reorganizing it wasn't worth the trouble.
		*/
	};

	// The size of the permutation
	int _size;

	// The divisors
	std::vector<int> _divisors;

	// The number of existing divisors 
	int _numDivisors;

	// The maximum number of cycles
	int _maxCycles;

	// The current number of cycles
	int _numCycles;

	// The indices used
	std::vector<bool> _used;

	// The number of indices participating in cycles
	int _numUsed;

	// The current cycle structure of the permutation
	std::vector<Cycle*> _cycles;

	// Is this the first permutation
	bool _firstPermutation;

	// The current permutation
	std::vector<int> _index;

	// operation order
	int _operationOrder;

	/**
	* Locate the next free index starting at a given index, or -1 if not found
	*/
	int findFreeIndex(int startingAt);

	/**
	* Find the divisors of the operation order and put them into the array
	*/
	void findDivisors();

	/**
	* Cancel a cycle - return all indices participating in it to
	* the identity.
	*/
	void cancelCycle(Cycle *c, bool clearUsed);

	/**
	* This function simply applies the given cycle -
	* Causes the permutation to contain it
	*/
	void applyCycle(Cycle *c);

	int initCycle(Cycle* c, int divisorIndex, int startingIndex, int ignoreUsed);
	void freeCycle(Cycle *c);

	/**
	* Clear the used array
	*/
	void clearPerm();

	/**
	* Try to advance a cycle, return FALSE upon failure
	*/
	int advanceCycle(int cycleIndex);

public:
	/**
	* Initialize a new permuter
	*
	* @param size The size of the permutations to generate (1 .. size)
	* @param groupSize The size of allowed cycles (in addition to cycles of size 1)
	* @param addGroupsOfTwo Should cycles of size 2 also be allowed
	*
	* @return The new permuter object
	*/
	Permuter(int size, int groupSize, int addGroupsOfTwo);
	~Permuter();

	/**
	* Generate the next permutation
	*
	* @param p The permuter object
	*
	* @return TRUE if we have not yet reached the end of the enumerator, FALSE otherwise
	*/
	bool next();

	/**
	* Reset the permuter to re-initialize enumeration
	*/
	void reset();

	/*
	* Returns one permuation element
	*
	* @param index The element index in the permutation
	* @returns an element at index @index
	*/
	const int elementAt(int index) const
	{
		return (*this)[index];
	}

	const int operator[](int index) const;
};
#endif
