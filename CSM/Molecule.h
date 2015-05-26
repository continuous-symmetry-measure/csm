/*
 * Author: shadi lahham
 *
 * A simple molecule structure with :
 *  - atom positions XYZ
 *  - atom symbols H|C|Na|Cl .. etc
 *  - adjacency matrix to represent connectivity 1|0
 *  - valency of each atom
 *  - similarity: d-similar atoms have the same symbol and same structure
 *                of children and grandchildren up to a depth of d
 *
 * Converted to C++ by Itay Zandbank
 */

#ifndef MOLECULE_H
#define MOLECULE_H

#include "csmlib.h"
#include <openbabel/mol.h>

#ifndef SQR
#define SQR(x)      ((x) * (x))
#endif

class Molecule
{
private:
	int _size;
	double** _pos;           // atom positions XYZ
	char** _symbol;          // atom symbols
	int** _adjacent;         // represent connectivity
	std::vector<int> _valency;           // valency of each atom
	std::vector<int> _similar;           // similarity
	int  _groupNum;          // the number of groups of similarity
	double _norm;	     // The normalization factor
	std::vector<double> _mass;	     // The atomic masses

private:
	Molecule(size_t size);  // Private constructor forces creation through the factory methods

	void replaceSymbols();
	void initSimilarity(int depth);
	int isSimilar(int a, int b);
	Molecule* copy(int* selectedAtoms, int selectedAtomsSize, bool updateSimilarity);

public:
	~Molecule();
	static Molecule *create(FILE *in, FILE *err, bool replaceSym);
	static Molecule* createFromOBMol(OpenBabel::OBMol &obmol, bool replaceSym, bool useMass = false);
	static Molecule* createFromPython(const std::vector<python_atom> &atoms);

public:
	int getGroup(int num, int* buff);
	int getGroupSize(int num);
	int getMaxGroupSize();
	Molecule* stripAtoms(char** removeList, int removeListSize, int updateSimilarity);
	bool normalizeMolecule(bool keepCenter);

	void print();
	void printBasic();
	void printSimilar();
	void printDebug();
	void printDebug2();

	/* Get properties */
	int groupNum() const { return _groupNum; }
	int size() const { return _size; }
	char *symbol(int index) const { return _symbol[index]; }
	char **symbols() const { return _symbol;  }
	double norm() const { return _norm; }
	double mass(int index) const { return _mass[index]; }
	int similar(int index) const { return _similar[index]; }
	int adjacent(int i, int j) const { return _adjacent[i][j]; }
	int valency(int index) const { return _valency[index]; }

	// _pos is exposed as is, because too much code expects a 2D matrix here, and changing it
	// would require too much modifications.
	// TODO: Rethink this once we replace nrutil
	double **pos() const { return _pos; }
};

#endif
