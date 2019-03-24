/*
 * Implementation of the math wrapper routines
 *
 * Written by Itay Zandbank
 */

#include "math_wrappers.h"
#include <iostream>
#include <vector>
#include <complex>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#if EIGEN_WORLD_VERSION==3 && EIGEN_MAJOR_VERSION==3 && EIGEN_MINOR_VERSION < 3
#error Unsupported Eigen Version, please upgrade to Eigen 3.3.3 or newer
#endif
// from rpoly.c Jenkins-Traub Polynomial Solver
extern "C" {
	int rpoly(double *op, int degree, double *zeror, double *zeroi);
}

/*
 * Find the roots of a polynomial specified by coefficients.
 * coefficients[0] is the coefficient of the highest degree, coefficients[1] of the second
 * highest and so on.
 *
 * This function delegates the calculations to rpoly.c
 */
std::vector<std::complex<double> > FindPolyRoots(const std::vector<double>& coefficients)
{
	// Prepare all the rpoly arguments
	// Allocate all the necessary memory manually, and copy the coefficients, since rpoly
	// might change them - who knows.
	size_t highest_degree = coefficients.size() - 1;
	double *zeror = new double[highest_degree];  // Not using std::vector<double>.data() because it just seems wrong
	double *zeroi = new double[highest_degree];  // although it is quite acceptable http://stackoverflow.com/questions/18759692/stdvector-write-directly-to-the-internal-array
	double *coeffs = new double[highest_degree + 1];

	for (int i = 0; i < highest_degree + 1; i++)
		coeffs[i] = coefficients[i];

	rpoly(coeffs, highest_degree, zeror, zeroi);
		
	delete[] coeffs;

	std::vector<std::complex<double> > roots;
	for (int i = 0; i < highest_degree; i++)
		roots.push_back(std::complex<double>(zeror[i], zeroi[i]));

	delete[] zeror;
	delete[] zeroi;

	return roots;
}


/*
 * Return the eigenvectors and eigenvalues of a 3x3 matrix.
 *
 * CSM only uses 3x3 matrices, so there was no point in supporting other matrix sizes.
 *
 * This is a thin wrapper around Eigen's EigenSolver.
 */

std::vector<EigenResult> GetEigens(const double matrix[3][3])
{
	//this function is not called by the python code
	Eigen::Matrix3d m;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
		{
			m(i, j) = matrix[i][j];
		}

	Eigen::EigenSolver<Eigen::Matrix3d> solver(m, true);
	std::vector<EigenResult> results(3);
	for (int i = 0; i < 3; i++)
	{
		results[i].value = solver.eigenvalues()[i].real();
		results[i].vector.resize(3);
		for (int j = 0; j < 3; j++)
		{
			results[i].vector[j] = solver.eigenvectors().col(i)[j].real();
		}
	}

	return results;
}

void GetEigens(const double matrix[3][3], double eigenVectors[3][3], double eigenValues[3])
{
	Eigen::Matrix3d m;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
		{
			//std::cout << "matij: " << matrix[i][j] << "\n";
			m(j,i) = matrix[i][j];
			//std::cout << "mji: " << m(j,i) << "\n";
		}

	Eigen::EigenSolver<Eigen::Matrix3d> solver(m, true);

	for (int i = 0; i < 3; i++)
	{
		eigenValues[i] = solver.eigenvalues()[i].real();
		//std::cout << "eigenvalues i: " << eigenValues[i] << "\n";
		for (int j = 0; j < 3; j++)
		{
			eigenVectors[i][j] = solver.eigenvectors().col(i)[j].real();
		}
	}
}
