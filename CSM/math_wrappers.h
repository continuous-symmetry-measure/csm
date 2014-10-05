/*
 * Mathematical routines for CSM
 *
 * These are thin wrappers around other code, used to move from the old NR code to new open-source code.
 */

#ifndef MATH_WRAPPERS_H
#define MATH_WRAPPERS_H

#include <vector>
#include <complex>

/* Returns the root of the polynomial passed in coefficients. */
std::vector<std::complex<double>> FindPolyRoots(const std::vector<double>& coefficients);

struct EigenResult
{
	std::vector<double> vector;
	double value;
};
std::vector<EigenResult> GetEigens(const double matrix[3][3]);

#endif