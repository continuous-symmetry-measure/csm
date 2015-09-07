/*
 * The main CSM calculations.
 *
 * Extracted from the old mainRot.cpp during the 2014 organization.
 * Created by Itay Zandbank
 *
 */

#ifndef CALCULATIONS_H
#define CALCULATIONS_H

#include "options.h"
#include "Molecule.h"

#define MAXDOUBLE  100000000.0
#define MINDOUBLE  1e-8
#define GROUPSIZE_LIMIT 15
#define GROUPSIZE_FACTOR 1.32e11
#define APPROX_RUN_PER_SEC 8e4
#define ZERO_IM_PART_MAX (1e-3)
#define MIN_GROUPS_FOR_OUTLIERS 10

void runSinglePerm(Molecule* m, double** outAtoms, int* perm, double* csm, double* dir, double* dMin, OperationType type);
void findBestPerm(Molecule* m, double** outAtoms, int* optimalPerm, double* csm, double* dir, double* dMin, OperationType type);
void findBestPermUsingDir(Molecule* m, double** outAtoms, int* optimalPerm, double* csm, double* dir, double* dMin, OperationType type);
void findSymmetryDirection(Molecule *m, double  ***dirs, int *n_dirs, OperationType type);
void estimatePerm(Molecule* m, int *perm, double *dir, OperationType type);
void lineFit(double **points, int nPoints, double **dirs, int* outlies);
void planeFit(double **points, int nPoints, double **dirs, int *outliers);
void initIndexArrays(Molecule* m, int* posToIdx, int* idxToPos);
double computeLocalCSM(Molecule* m, double *localCSM, int *perm, double *dir, OperationType type);

double calcRefPlane(Molecule* m, int* perm, double *dir, OperationType type);
double createSymmetricStructure(Molecule* m, double **outAtoms, int *perm, double *dir, OperationType type, double dMin);

#endif