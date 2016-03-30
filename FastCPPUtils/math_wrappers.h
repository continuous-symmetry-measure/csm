/*
 * Mathematical routine wrappers for CSM
 *
 * CSM needs to find roots of polynomials and the eigenvectors and eigenvalues of a matrix.
 * This was originally done by various third-party routines that were used directly in code.
 * 
 * During transition to C++, a more C++ friendly interface was designed. In addition, the
 * underlying routines have been hidden behind these two simple functions. We've used this 
 * abstraction to move from closed-source Numerical Recipes code to the Open Source Eigen library.
 */

#ifndef MATH_WRAPPERS_H
#define MATH_WRAPPERS_H

#include <vector>
#include <complex>

/* Returns the root of the polynomial passed in coefficients. */
std::vector<std::complex<double> > FindPolyRoots(const std::vector<double>& coefficients);

/*
 * An eigenvalue, eigenvector pair.
 */
struct EigenResult
{
	double value;
	std::vector<double> vector;
};

/*
 * Returns all the eigenvalues and eigenvectors of a 3x3 matrix.
 */
std::vector<EigenResult> GetEigens(const double matrix[3][3]);

#endif