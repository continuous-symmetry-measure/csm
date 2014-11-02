/*
 * Some utility mathematical functions.
 *
 * Extracted from mainRot.cpp during the 2014 reorganization.
 * Created by Itay Zandbank
 */

#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include "Molecule.h"

void normalize(double **coords, Molecule *m);
double Magnitude(double *Point1, double *Point2);
double computeDistanceFromLine(double *Point, double *LineStart, double *LineEnd);
double findMedian(double arr[], int n);


#endif