#ifndef BABEL_ADAPTER_H
#define BABEL_ADAPTER_H

// Include Open Babel classes for OBMol and OBConversion
extern "C" {
#include <openbabel/mol.h>
}
#include "Molecule.h"

using namespace OpenBabel;

/** 
 * Get the atomic mass of a given element
 * @param atomName The atom's symbol
 *  
 * @return The atomic mass
 */ 
double getAtomicMass(char *atomName);

/**
 * Updates the coordinates of the OpenBabel Molecule according to the Molecule data
 * 
 * @param obmol The OpenBable molecule
 * @param outAtoms The output atoms' coordinates
 */
void updateCoordinates(OBMol& obmol, double **outAtoms);

/**
 * Read a molecule from a file file
 * 
 * @param filename The file name to read from
 * @param format   The file format in Open Babel terms (if NULL is given, attempt to get format 
 * 	  from the suffix)
 * @param babelBond Should the convertor let openbabel choose bonding
 * @return The OBMol read from the file
 */
OBMol readMolecule (char *filename, const std::string format, int babelBond);

/**
 * Write a mol to file
 *  
 * @param mol The molecule
 * @param format The output format
 * @param file The file to write output to
 * @param filename The file's name (for extension-finding purpose)
 */
void writeMolecule(OBMol& mol, const std::string format, FILE* file, char *filename);

#endif
