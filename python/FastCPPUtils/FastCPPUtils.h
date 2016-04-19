/*
 * A few functions that are a lot faster in C++ than their Numpy equivalent
 * By Itay Zandbank
 */

#ifndef FASTCPPUTILS_H
#define FASTCPPUTILS_H

extern "C"
{
	int rpoly(double *op, int degree, double *zeror, double *zeroi);
}
void GetEigens(const double matrix[3][3], double eigenVectors[3][3], double eigenValues[3]);

#endif
