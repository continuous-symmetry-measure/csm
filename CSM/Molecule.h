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
 */

#ifndef MOLECULE_H
#define MOLECULE_H

#include <openbabel/mol.h>

#ifndef SQR
#define SQR(x)      ((x) * (x))
#endif

class Molecule
{
public:
    int _size;
    double** _pos;           // atom positions XYZ
    char** _symbol;          // atom symbols
    int** _adjacent;         // represent connectivity
    int* _valency;           // valency of each atom
    int* _similar;           // similarity
    int* _marked;            // for marking atoms - general use
    int  _groupNum;          // the number of groups of similarity
    double _norm;	     // The normalization factor
    double* _mass;	     // The atomic masses

private:
	Molecule(int size);  // Private constructor forces creation through the factory methods

	void replaceSymbols();
	void initSimilarity(int depth);
	void setMarked(int state);
	int isSimilar(int a, int b);

public:
	~Molecule();
	static Molecule *create(FILE *in, FILE *err, bool replaceSym);
	static Molecule *createPDB(FILE *in, FILE *err, bool replaceSym);
	static Molecule* createFromOBMol(OpenBabel::OBMol &obmol, bool replaceSym, bool useMass = false);
	Molecule *copy(int *selectedAtoms, int selectedAtomsSize, bool updateSimilarity);

	void print();
	void printBasic();
	void printSimilar();
	void printDebug();
	void printDebug2();


};

int getGroup(Molecule *m,int num,int* buff);

int getGroupSize(Molecule *m,int num);

int getMaxGroupSize(Molecule *m);

Molecule* stripAtoms(Molecule *m, char** removeList, int removeListSize, int updateSimilarity);

int normalizeMolecule(Molecule *m, bool keepCenter);
#endif
