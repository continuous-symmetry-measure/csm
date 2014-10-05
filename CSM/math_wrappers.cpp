/*
 * Implementation of the math wrapper routines
 */

#include "math_wrappers.h"

#include <vector>
#include <complex>
#include "math_utils.h"
#include "logging.h"
#include "dmatrix.h"
#include "dvector.h"

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

std::vector<EigenResult> GetEigens(const double matrix[3][3])
{
	csm_utils::dmatrix tmatrix(1, 3, 1, 3);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			tmatrix[i+1][j+1] = matrix[i][j];
	csm_utils::dvector diag(1, 3), secdiag(1, 3);

	tred2(tmatrix, 3, diag, secdiag);
	tqli(diag, secdiag, 3, tmatrix);

	std::vector<EigenResult> results(3);
	for (int i = 0; i < 3; i++)
	{
		results[i].value = diag[i+1];
		results[i].vector.resize(3);
		for (int j = 0; j < 3; j++)
			results[i].vector[j] = tmatrix[j+1][i+1];
	}

	return results;
}