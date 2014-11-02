/*
 * Various print-out utility functions.
 *
 * Taken from mainRot.cpp during the 2014 organization.
 * Created by Itay Zandbank
 */

#ifndef PRINTOUTS_H
#define PRINTOUTS_H

#include "Molecule.h"
#include "babelAdapter.h"
#include <string>

extern char opName[100];
extern bool printNorm;
extern bool printLocal;
extern bool writeOpenu;
extern std::string format;

void printOutput(Molecule* m, double** outAtoms, double csm, double *dir, double dMin, FILE *out, double* localCSM);
void printPDBATOM(Molecule* m, FILE* f, char** sym, double** pos);
void printOutputPDB(Molecule* m, double** outAtoms, double csm, double *dir, double dMin, FILE *out);
void printOutputFormat(Molecule* m, OBMol& mol, double** outAtoms, double csm, double *dir, double dMin, FILE *out, char *fname, double* localCSM);

#endif