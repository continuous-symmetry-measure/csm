#include <iostream>
#include <sstream>
#include <string.h>
#include <limits.h>

// Include Open Babel classes for OBMol and OBConversion
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/obconversion.h>
#include <stdio.h>
#include "babelAdapter.h"

using namespace OpenBabel;
using namespace std;

#include "Molecule.h"
#include "elements.h"

extern void initSimilarity(Molecule *m,int depth);		
extern void replaceSymbols(Molecule* m);

/** 
 * Get the atomic mass of a given element
 * @param atomName The atom's symbol
 *  
 * @return The atomic mass
 */ 
double getAtomicMass(char *atomName) { 
	for (int i = 0; i < ELEMENTS.size(); i++) { 
		if (!strcmp(atomName, ELEMENTS[i].c_str())) { 
			OBAtom a;
			a.SetAtomicNum(i+1);
			a.SetIsotope(0);
			return a.GetAtomicMass();
		}
	}
	printf("Failed: Unknown atom type %s\n", atomName);
	exit(1);
	return 1;
} 

/**
 * Updates the coordinates of the OpenBabel Molecule according to the Molecule data
 * 
 * @param obmol The OpenBable molecule
 * @param outAtoms The output atoms' coordinates
 */
void updateCoordinates(OBMol& obmol, double **outAtoms) {
	int numAtoms = obmol.NumAtoms();
	int i;
	for (i = 0; i < numAtoms; i++) {
		OBAtom* atom = obmol.GetAtom(i + 1);
		atom->SetVector(outAtoms[i][0], outAtoms[i][1], outAtoms[i][2]);
	}
}

/**
 * Read a molecule from a file file
 * 
 * @param filename The file name to read from
 * @param format   The file format in Open Babel terms (if NULL is given, attempt to get format 
 * 	  from the suffix)
 * @param babelBond Should the convertor let openbabel choose bonding
 * @return The OBMol read from the file
 */
OBMol readMolecule (char *filename, const char *format, int babelBond) {
	OBMol mol;
	OBConversion conv;		
	if (format == NULL) {
		OBFormat* f = OBConversion::FormatFromExt(filename);
		if (f == NULL) {	
			printf("Error discovering format from filename %s\n", filename);
			exit(1);
		}			
		if (!conv.SetInFormat(f)) {
			printf ("Error setting openbabel format\n");
			exit(1);
		}
	} else {
		if (!conv.SetInFormat(format)) {
			printf ("Error setting input format to %s\n", format);
			exit(1);
		}
	}
	if (!babelBond) conv.SetOptions("b", OBConversion::INOPTIONS);
	if (!conv.ReadFile(&mol, filename)) {
		printf ("Error reading file %s using OpenBabel\n", filename);
		exit(1);		
	}		
	return mol;
}

/**
 * Write a mol to file
 *  
 * @param mol The molecule
 * @param format The output format
 * @param file The file to write output to
 * @param filename The file's name (for extension-finding purpose)
 */
void writeMolecule(OBMol& mol, const char *format, FILE* file, char *filename) {
	ostringstream os;
	OBConversion conv;
	if (format == NULL) {
		OBFormat* f = OBConversion::FormatFromExt(filename);
		if (f != NULL) {	
			printf("Error discovering format from filename %s\n", filename);
			exit(1);
		}
		if (!conv.SetOutFormat(f)) {
			printf ("Error setting openbabel format\n");
			exit(1);
		}		
	} else {
		if (!conv.SetOutFormat(format)) {
			printf ("Error setting output format to %s\n", format);
			exit(1);
		}
	}
	if (!conv.Write(&mol, &os)) {
		printf ("Error writing data file using OpenBabel\n");
		exit(1);
	}
	fprintf(file,"%s\n",os.str().c_str());
}
