/*
 * Author: shadi lahham
 *
 * contains helper functions for the C2 operation
 *
 */

#ifndef C2HELPERS_H
#define C2HELPERS_H

#include "Molecule.h"

#define MAXDOUBLE  100000000.0
#define MINDOUBLE  1e-8
#define PI       3.14159265358979323846


double** createAntimer(Molecule* m);
void freeAntimer(Molecule* m,double** antimer);
double chiralityFunction(double *parray[], double *anti_parray[],
                        int num_of_atoms, double *dir, double *d_min, double eulerParam[]);
void foldingUnfolding(double *parray[], double *anti_parray[],double *uarray[],
                      int num_of_atoms, double *dir_cos,double d_min, double eulerParam[]);


#endif
