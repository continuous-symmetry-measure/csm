/*
 * Implementation of the math wrapper routines
 */

#include "math_wrappers.h"

#include <vector>
#include <complex>

// from rpoly.c Jenkins-Traub Polynomial Solver
extern "C" {
	int rpoly(double *op, int degree, double *zeror, double *zeroi);
}

std::vector<std::complex<double>> FindPolyRoots(const std::vector<double>& coefficients)
{
	// First version, using rpoly
	int highest_degree = coefficients.size() - 1;
	double *zeror = new double[highest_degree];
	double *zeroi = new double[highest_degree];
	double *coeffs = new double[highest_degree + 1];

	for (int i = 0; i < highest_degree + 1; i++)
		coeffs[i] = coefficients[i];

	rpoly(coeffs, highest_degree, zeror, zeroi);
		
	delete[] coeffs;

	std::vector<std::complex<double>> roots;
	for (int i = 0; i < highest_degree; i++)
		roots.push_back(std::complex<double>(zeror[i], zeroi[i]));

	delete[] zeror;
	delete[] zeroi;

	return roots;
}