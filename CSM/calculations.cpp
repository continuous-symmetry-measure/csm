/*
* The main CSM calculations.
*
* Extracted from the old mainRot.cpp during the 2014 organization.
* Created by Itay Zandbank
*
*/

#include <boost/format.hpp>
#include <sstream>
#include <iomanip>
#include <complex>
#include <memory.h>

#include "calculations.h"
#include "dmatrix.h"
#include "dvector.h"
#include "options.h"
#include "math_wrappers.h"
#include "groupPermuter.h"
#include "permuter.h"
#include "math_utils.h"
#include "logging.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

struct distRecord {
	double distance;
	int row;
	int col;
};

int distComp(const void *pd1, const void* pd2) {
	const struct distRecord* d1 = (const struct distRecord*)pd1;
	const struct distRecord* d2 = (const struct distRecord*)pd2;
	if (d1->distance > d2->distance) {
		return 1;
	}
	else if (d1->distance < d2->distance) {
		return -1;
	}
	return 0;
}


/*
* creates position to index , and index to position translation arrays
*/
void initIndexArrays(Molecule* m, int* posToIdx, int* idxToPos){

	int i, j, counter;

	counter = 0;

	// build idxToPos
	for (j = 1; j <= m->groupNum(); j++){
		for (i = 0; i< m->size(); i++){
			if (m->similar(i) == j){
				idxToPos[counter] = i;
				counter++;
			}
		}
	}

	// build posToIdx
	for (i = 0; i< m->size(); i++){
		posToIdx[idxToPos[i]] = i;
	}

}

/**
* Compute the part of the A Matrix relevant to the current permutation
*/
void computeMatrix(double **parray, int *perm, int size, double(*coef)[3][3], double multiplier) {
	int atom;
	for (atom = 0; atom < size; ++atom) {
		(*coef)[0][0] += 2.0 * parray[atom][0] * parray[perm[atom]][0] * multiplier;
		(*coef)[1][1] += 2.0 * parray[atom][1] * parray[perm[atom]][1] * multiplier;
		(*coef)[2][2] += 2.0 * parray[atom][2] * parray[perm[atom]][2] * multiplier;
		(*coef)[1][0] += (parray[atom][0] * parray[perm[atom]][1] +
			parray[atom][1] * parray[perm[atom]][0]) * multiplier;
		(*coef)[2][0] += (parray[atom][0] * parray[perm[atom]][2] +
			parray[atom][2] * parray[perm[atom]][0]) * multiplier;
		(*coef)[2][1] += (parray[atom][1] * parray[perm[atom]][2] +
			parray[atom][2] * parray[perm[atom]][1]) * multiplier;
	}
	(*coef)[0][1] = (*coef)[1][0];
	(*coef)[0][2] = (*coef)[2][0];
	(*coef)[1][2] = (*coef)[2][1];
}


/**
* Compute the part of the B vec relevant to the current permutation
*/
void computeVector(double **parray, int *perm, int size, double(*vec)[3], double multiplier) {
	int atom;
	for (atom = 0; atom < size; ++atom) {
		(*vec)[0] += multiplier * (parray[atom][1] * parray[perm[atom]][2] - parray[atom][2] * parray[perm[atom]][1]);
		(*vec)[1] += multiplier * (parray[atom][2] * parray[perm[atom]][0] - parray[atom][0] * parray[perm[atom]][2]);
		(*vec)[2] += multiplier * (parray[atom][0] * parray[perm[atom]][1] - parray[atom][1] * parray[perm[atom]][0]);
	}
}

/**
* Calculate the best axis, and compute the CSM for it, given a pairing of the indices (perm)
*/
double calcRefPlane(Molecule* m, int* perm, double *dir, OperationType type) {
	csm_utils::dmatrix copyMat(1, 3, 1, 3);
	csm_utils::dvector diag(1, 3);
	csm_utils::dvector temp(1, 3);

	double matrix[3][3] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };
	double vec[3] = { 0.0, 0.0, 0.0 };
	double scalar[3];
	double maxval, scl, angle;
	std::vector<double> coeffs(7);
	int i, j;
	int *curPerm = (int *)malloc(sizeof(int) * m->size());
	double csm, dists;

	int isImproper = (type != CN) ? true : false;
	int isZeroAngle = (type == CS) ? true : false;

	LOG(debug) << "calcRefPlane called";
	stringstream permstrm;
	for (int i = 0; i < m->size(); i++)
	{
		if (i > 0)
			permstrm << ", ";
		permstrm << perm[i];
	}
	LOG(debug) << "Permutation is " << permstrm.str();
	LOG(debug) << "Direction is " << setprecision(2) << fixed << dir[0] << " " << dir[1] << " " << dir[2];

	// initialize identity permutation
	for (i = 0; i < m->size(); i++) {
		curPerm[i] = i;
	}

	// compute matrices according to current perm and its powers (the identity does not contribute anyway)
	for (i = 1; i < options.opOrder; i++) {
		angle = isZeroAngle ? 0.0 : (2 * M_PI * i / options.opOrder);
		// The i'th power of the permutation
		for (j = 0; j < m->size(); j++) {
			curPerm[j] = perm[curPerm[j]];
		}
		if (isImproper && ((i % 2) == 1)) {
			computeMatrix(m->pos(), curPerm, m->size(), &matrix, -1 - cos(angle));
		}
		else {
			computeMatrix(m->pos(), curPerm, m->size(), &matrix, 1 - cos(angle));
		}
		computeVector(m->pos(), curPerm, m->size(), &vec, sin(angle));
	}

	LOG(debug).unsetf(ios_base::fixed);
	LOG(debug) << "Computed matrix is:" << setprecision(4);
	LOG(debug) << matrix[0][0] << " " << matrix[0][1] << " " << matrix[0][2];
	LOG(debug) << matrix[1][0] << " " << matrix[1][1] << " " << matrix[1][2];
	LOG(debug) << matrix[2][0] << " " << matrix[2][1] << " " << matrix[2][2];

	vector<EigenResult> eigens = GetEigens(matrix);

	// compute square of scalar multiplications of eigen vectors with b
	for (i = 0; i < 3; i++) {
		scalar[i] = 0.0;
		for (j = 0; j < 3; j++) {
			scalar[i] += vec[j] * eigens[i].vector[j]; // copyVec[j + 1] * copyMat[j+1][i+1];
			copyMat[j + 1][i + 1] = eigens[i].vector[j];
		}
		temp[i + 1] = scalar[i] * scalar[i];
	}
	for (i = 0; i < 3; i++)
		diag[i + 1] = eigens[i].value;

	// build the polynomial
	coeffs[0] = 1.0;	// x^6
	coeffs[1] = -2 * (diag[1] + diag[2] + diag[3]);	// x^5
	coeffs[2] = diag[1] * diag[1] + diag[2] * diag[2] + diag[3] * diag[3] -
		temp[1] - temp[2] - temp[3] +
		4 * (diag[1] * diag[2] + diag[1] * diag[3] + diag[2] * diag[3]); // x^4
	coeffs[3] = -8 * diag[1] * diag[2] * diag[3] +
		2 * (temp[1] * diag[2] +
		temp[1] * diag[3] +
		temp[2] * diag[1] +
		temp[2] * diag[3] +
		temp[3] * diag[1] +
		temp[3] * diag[2] -
		diag[1] * diag[3] * diag[3] -
		diag[1] * diag[1] * diag[2] -
		diag[1] * diag[1] * diag[3] -
		diag[1] * diag[2] * diag[2] -
		diag[2] * diag[2] * diag[3] -
		diag[2] * diag[3] * diag[3]); // x^3
	coeffs[4] = 4 *
		((diag[1] * diag[2] * diag[3] * (diag[1] + diag[2] + diag[3]) -
		(temp[3] * diag[1] * diag[2] +
		temp[2] * diag[1] * diag[3] +
		temp[1] * diag[3] * diag[2]))) -
		temp[1] * (diag[2] * diag[2] + diag[3] * diag[3]) -
		temp[2] * (diag[1] * diag[1] + diag[3] * diag[3]) -
		temp[3] * (diag[1] * diag[1] + diag[2] * diag[2]) +
		diag[1] * diag[1] * diag[2] * diag[2] +
		diag[2] * diag[2] * diag[3] * diag[3] +
		diag[1] * diag[1] * diag[3] * diag[3]; // x^2
	coeffs[5] = 2 *
		(temp[1] * diag[2] * diag[3] * (diag[2] + diag[3]) +
		temp[2] * diag[1] * diag[3] * (diag[1] + diag[3]) +
		temp[3] * diag[1] * diag[2] * (diag[1] + diag[2]))
		- 2 *
		(diag[1] * diag[2] * diag[2] * diag[3] * diag[3] +
		diag[1] * diag[1] * diag[2] * diag[3] * diag[3] +
		diag[1] * diag[1] * diag[2] * diag[2] * diag[3]); // x
	coeffs[6] = -temp[1] * diag[2] * diag[2] * diag[3] * diag[3] -
		temp[2] * diag[1] * diag[1] * diag[3] * diag[3] -
		temp[3] * diag[1] * diag[1] * diag[2] * diag[2] +
		diag[1] * diag[1] * diag[2] * diag[2] * diag[3] * diag[3]; // 1

	// solve polynomial and find maximum eigenvalue and eigen vec
	LOG(debug) << "Coefficients: " << coeffs[0] << ", " << coeffs[1] << ", " << coeffs[2] << ", " << coeffs[3] << ", " << coeffs[4] << ", " << coeffs[5] << ", " << coeffs[6];
	vector<complex<double> > roots = FindPolyRoots(coeffs);
	// rpoly(coeffs, 6, rtr, rti);

	LOG(debug) << "rtr: " << roots[0] << " " << roots[1] << " " << roots[2] << " " << roots[3] << " " << roots[4] << " " << roots[5];

	maxval = -MAXDOUBLE;
	for (i = 0; i < 6; i++) {
		if (maxval < roots[i].real() && (fabs(roots[i].imag()) < ZERO_IM_PART_MAX))
			maxval = roots[i].real();
	}

	LOG(debug) << setprecision(6) << fixed << "diag: " << diag[1] << " " << diag[2] << " " << diag[3];
	scl = 0.0;
	if ((isZeroAngle) || (options.opOrder == 2)) {
		// If we are in zero angle case, we should pick the direction matching maxval
		double minDist = MAXDOUBLE;
		int minarg = 0;
		for (i = 1; i <= 3; i++) {
			if (fabs(diag[i] - maxval) < minDist) {
				minDist = fabs(diag[i] - maxval);
				minarg = i;
			}
		}
		for (i = 0; i < 3; i++) {
			dir[i] = copyMat[i + 1][minarg];
		}
	}
	else {
		for (i = 1; i <= 3; i++) {
			dir[i - 1] = 0.0;
			for (j = 1; j <= 3; j++) {
				// error safety
				if (fabs(diag[j] - maxval) < 1e-6) {
					dir[i - 1] = copyMat[i][j];
					break;
				}
				else {
					dir[i - 1] += scalar[j - 1] / (diag[j] - maxval) * copyMat[i][j];
				}
			}
			scl += dir[i - 1] * vec[i - 1];
		}
	}

	// initialize identity permutation
	for (i = 0; i < m->size(); i++) {
		curPerm[i] = i;
	}

	// compute CSM	
	csm = 0.0;
	for (i = 0; i < options.opOrder; i++) {
		// This can be more efficient - save results of matrix computation
		if (i != 0) {
			angle = isZeroAngle ? 0.0 : ((2 * M_PI * i) / options.opOrder);
			dists = 0.0;
			for (j = 0; j < m->size(); j++) {
				// i'th power of permutation
				curPerm[j] = perm[curPerm[j]];
			}
			for (j = 0; j < m->size(); j++) {
				dists += (m->pos()[j][0] * m->pos()[curPerm[j]][0] +
					m->pos()[j][1] * m->pos()[curPerm[j]][1] +
					m->pos()[j][2] * m->pos()[curPerm[j]][2]);
			}
			csm += cos(angle) * dists;

		}
		else {
			csm += 1.0;
		}
	}


	LOG(debug) << setprecision(6) << fixed << "csm=" << csm << " maxval=" << maxval << " scl=" << scl;
	LOG(debug) << setprecision(6) << fixed << "dir: " << dir[0] << " " << dir[1] << " " << dir[2];
	csm += (maxval - scl) / 2;
	csm = fabs(100 * (1.0 - csm / options.opOrder));
	free(curPerm);

	LOG(debug) << setprecision(6) << fixed << "dir - csm: " << dir[0] << " " << dir[1] << " " << dir[2] << " - " << csm;

	return csm;
}

double createSymmetricStructure(Molecule* m, double **outAtoms, int *perm, double *dir, OperationType type, double dMin) {
	int isImproper = (type != CN) ? true : false;
	int isZeroAngle = (type == CS) ? true : false;
	int i, j, k, l;
	int *curPerm = (int *)malloc(sizeof(int) * m->size());
	double rotaionMatrix[3][3];
	double tmpMatrix[3][3] = { { 0.0, -dir[2], dir[1] }, { dir[2], 0.0, -dir[0] }, { -dir[1], dir[0], 0.0 } };
	double angle;
	double res = 0.0;

	for (i = 0; i < m->size(); i++) {
		// initialize with identity operation
		curPerm[i] = i;
		for (j = 0; j < 3; j++) {
			outAtoms[i][j] = m->pos()[i][j];
		}
	}

	for (i = 1; i < options.opOrder; i++) {
		angle = isZeroAngle ? 0.0 : (2 * M_PI * i / options.opOrder);
		int factor = ((isImproper && (i % 2) == 1) ? (-1) : 1);
		for (j = 0; j < m->size(); j++) {	
			curPerm[j] = perm[curPerm[j]];
		}
		for (j = 0; j < 3; j++) {
			for (k = 0; k < 3; k++) {
				rotaionMatrix[j][k] =
					((j == k) ? cos(angle) : 0) +
					(factor - cos(angle)) * dir[j] * dir[k] +
					sin(angle) * tmpMatrix[j][k];
			}
		}

		for (j = 0; j < m->size(); j++) {
			for (k = 0; k < 3; k++) {
				for (l = 0; l < 3; l++) {
					outAtoms[j][k] += rotaionMatrix[k][l] * m->pos()[curPerm[j]][l];
					
				}
			}
		}
	}

	for (j = 0; j < m->size(); j++) {
		for (k = 0; k < 3; k++) {
			outAtoms[j][k] /= options.opOrder;
			outAtoms[j][k] *= dMin;
			res += SQR(outAtoms[j][k]);
			
		}
	}

	free(curPerm);

	return sqrt(res);
}

double computeLocalCSM(Molecule* m, double *localCSM, int *perm, double *dir, OperationType type) {
	int isImproper = (type != CN) ? true : false;
	int isZeroAngle = (type == CS) ? true : false;
	int i, j, k, l;
	int *curPerm = (int *)malloc(sizeof(int) * m->size());
	double rotaionMatrix[3][3];
	double tmpMatrix[3][3] = { { 0.0, -dir[2], dir[1] }, { dir[2], 0.0, -dir[0] }, { -dir[1], dir[0], 0.0 } };
	double angle;

	double tempCSM = 0.0;


	for (i = 0; i < m->size(); i++) {
		// initialize with identity operation
		curPerm[i] = i;
		localCSM[i] = 0;
	}

	for (i = 1; i < options.opOrder; i++) {
		angle = isZeroAngle ? 0.0 : (2 * M_PI * i / options.opOrder);
		int factor = ((isImproper && (i % 2) == 1) ? (-1) : 1);
		for (j = 0; j < m->size(); j++) {
			curPerm[j] = perm[curPerm[j]];
		}
		for (j = 0; j < 3; j++) {
			for (k = 0; k < 3; k++) {
				rotaionMatrix[j][k] =
					((j == k) ? cos(angle) : 0) +
					(factor - cos(angle)) * dir[j] * dir[k] +
					sin(angle) * tmpMatrix[j][k];
			}
		}

		for (j = 0; j < m->size(); j++) {
			double rotated[3] = { 0.0, 0.0, 0.0 };
			for (k = 0; k < 3; k++) {
				for (l = 0; l < 3; l++) {
					rotated[k] += rotaionMatrix[k][l] * m->pos()[curPerm[j]][l];
				}
			}
			for (k = 0; k < 3; k++) {
				tempCSM += SQR(rotated[k] - m->pos()[j][k]);
				localCSM[j] += SQR(rotated[k] - m->pos()[j][k]);
			}
		}

	}

	for (j = 0; j < m->size(); j++) {
		localCSM[j] *= 100 / (2 * options.opOrder);
	}
	tempCSM *= 100.0 / (2 * options.opOrder);

	return tempCSM;
}

/*
* Calculates csm, dMin and directional cosines for a given permutation
*/
void runSinglePerm(Molecule* m, double** outAtoms, int *perm, double* csm, double* dir, double* dMin, OperationType type){
	*csm = calcRefPlane(m, perm, dir, type);

	// which is DMIN?
	*dMin = (1.0 - (*csm / 100 * options.opOrder / (options.opOrder - 1)));
	createSymmetricStructure(m, outAtoms, perm, dir, type, *dMin);
}

void findBestPermUsingDir(Molecule* m, double** outAtoms, int* perm, double* csm, double* dir, double* dMin, OperationType type) {
	estimatePerm(m, perm, dir, type);
	runSinglePerm(m, outAtoms, perm, csm, dir, dMin, type);
}

/**
* Finds an approximate permutation which can be used in the analytical computation.
*/
void findBestPerm(Molecule* m, double** outAtoms, int* perm, double* csm, double* dir, double* dMin, OperationType type) {
	LOG(debug) << "findBestPerm called";

	int *temp = (int*)malloc(sizeof(int) * m->size());
	int i = 0;
	// The algorithm aims to find the best perm which can then be used for the analytic solution	
	if ((type == CI) || (type == SN && options.opOrder == 2)) {
		// For inversion - simply find for each orbit the best matching - the closest after the operation.	
		// Do nothing - no need to find axis. yay.
		dir[0] = 1.0; dir[1] = 0.0; dir[2] = 0.0;
		estimatePerm(m, perm, dir, type);
		runSinglePerm(m, outAtoms, perm, csm, dir, dMin, type);
	}
	else {
		double** dirs;
		int n_dirs;
		int *bestPerm = (int*)malloc(sizeof(int) * m->size());
		double bestDir[3];
		int maxIters = 50;
		*csm = MAXDOUBLE;

		// Find an initial approximated symmetry axis/plain
		findSymmetryDirection(m, &dirs, &n_dirs, type);

		// Find the permutation best matching this direction - there are n_dirs results to the algorithm		
		for (i = 0; i < n_dirs; i++) {
			double dist = MAXDOUBLE, old = MAXDOUBLE, best = MAXDOUBLE;
			double tempDir[3];
			int iterNum = 1;
			estimatePerm(m, bestPerm, dirs[i], type);
			runSinglePerm(m, outAtoms, bestPerm, &dist, tempDir, dMin, type);
			memcpy(bestDir, tempDir, sizeof(double) * 3);
			old = MAXDOUBLE; best = dist;

			// solve analytically using this permutation, repeat until converged
			// Pick the best one
			// Stop once:
			// 1. The csm is less than 1e-4
			// 2. The difference between the old and new csm is less than 1%
			// 3. The max number of iterations has been reached
			while ((fabs(dist) > 1e-4) && (fabs(old - dist) / fabs(old) > 0.01) && (iterNum < maxIters)) {
				old = dist;
				estimatePerm(m, temp, tempDir, type);
				runSinglePerm(m, outAtoms, temp, &dist, tempDir, dMin, type);
				if (dist < best) {
					best = dist;
					memcpy(bestPerm, temp, sizeof(int) * m->size());
					memcpy(bestDir, tempDir, sizeof(double) * 3);
				}
				iterNum++;
				LOG(debug) << "Old csm: " << setprecision(4) << fixed << old << ", new csm " << dist;
			};

			// Keep the best solution so far...
			if (best < *csm) {
				*csm = best;
				memcpy(perm, bestPerm, sizeof(int) * m->size());
				memcpy(dir, bestDir, sizeof(double) * 3);
			}
			LOG(info) << "Attempt #" << i + 1 << ": best csm is " << setprecision(2) << fixed << best << " after " << iterNum << " iterations";
		}
		for (i = 0; i < n_dirs; i++) {
			free(dirs[i]);
		}
		free(dirs);
		free(bestPerm);

		// run once more to get the out atoms right !
		runSinglePerm(m, outAtoms, perm, csm, dir, dMin, type);
	}

#if 0 
	if (anneal) {

		// The initial probability is set to be 2 / N 
		// for a change of 0.0001
		double initialTemp = 0.0001 * log(1.0 * m->size() / 4.0);
		// printf("IT = %4.2f\n", initialTemp);
		double alpha = 0.9;
		int stepsPerPart = 2000;
		double t;
		double dist;
		double old = *csm;
		int factor = 1;
		double flipProb = 0.001;

		int **groups = (int**)malloc(sizeof(int *) * m->groupNum());
		int *groupSizes = (int*)malloc(sizeof(int) * m->groupNum());
		for (i = 0; i < m->groupNum(); i++) {
			groupSizes[i] = getGroupSize(m, i + 1);
			groups[i] = (int*)malloc(sizeof(int) * groupSizes[i]);
			getGroup(m, i + 1, groups[i]);
		}

		srand(time(NULL));
#ifndef _WIN32
		srand48(time(NULL));
#endif
		// try to anneal the best permutation
		memcpy(temp, perm, sizeof(int) * m->size());
		for (t = initialTemp; t >= 0.0001*initialTemp; t *= alpha) {
			for (i = 0; i < stepsPerPart /* * m->size() */; i++) {
				// select a node at random
				int first = rand() % m->size();

				// select another from its similarity group
				int second = groups
					[m->similar(first) - 1]
				[rand() % groupSizes[m->similar(first) - 1]];

				int temp1 = first, temp2 = second, firstSize = 1, secondSize = 1;
				if (first == second) continue;

				// find the orbit size of the two nodes, as well as the atoms preceeding them	
				while (temp[temp1] != first) {
					firstSize++;
					temp1 = temp[temp1];
				}
				while (temp[temp2] != second) {
					secondSize++;
					temp2 = temp[temp2];
				}
				// Now, finally, operate
				if (firstSize == 1 && secondSize == 1) {
					// both are single-orbits - so do nothing
				}
				else if (firstSize == 1) {
					// only first is a single orbit
					temp[temp2] = first;
					temp[first] = temp[second];
					temp[second] = second;
				}
				else if (secondSize == 1) {
					// only second is a single orbit
					temp[temp1] = second;
					temp[second] = temp[first];
					temp[first] = first;
				}
				else {
					// both are not single orbits
					int t = temp[first];
					temp[first] = temp[second];
					temp[second] = t;
					t = temp[temp1];
					temp[temp1] = temp[temp2];
					temp[temp2] = t;
				}

				runSinglePerm(m, outAtoms, temp, &dist, dir, dMin, type);
				if (factor == 1) {
					if (drand48() < flipProb) {
						factor = -1;
					}
				}
				else {
					// spend one tenth of the time in reverted state
					if (drand48() < (flipProb * 100)) {
						factor = 1;
					}
				}
				// printf("Tried to change from %4.2f to %4.2f: ", factor * old, factor * dist);
				if (dist < old || drand48() < exp(factor * (old - dist) / t)) {
					// if this improves - 

					/*if (dist < old) {
					printf("Better\n");
					} else {
					printf("Kept\n");
					} */
					if (dist < *csm) {
						*csm = dist;
						printf("Changed to %4.2f\n", *csm);
						memcpy(perm, temp, sizeof(int) * m->size());
					}
					old = dist;
				}
				else {
					// undo			
					//printf("Reverted\n");								
					if (firstSize == 1 && secondSize == 1) {
						// both are single-orbits - so do nothing
					}
					else if (firstSize == 1) {
						// only first is a single orbit
						temp[temp2] = second;
						temp[second] = temp[first];
						temp[first] = first;
					}
					else if (secondSize == 1) {
						// only second is a single orbit
						temp[temp1] = first;
						temp[first] = temp[second];
						temp[second] = second;
					}
					else {
						// both are not single orbits
						int t = temp[first];
						temp[first] = temp[second];
						temp[second] = t;
						t = temp[temp1];
						temp[temp1] = temp[temp2];
						temp[temp2] = t;
					}

				}
			}

			// printf("At T = %4.2f - csm %4.2f\n", t, *csm);
		}

		for (i = 0; i < m->groupNum(); i++) {
			free(groups[i]);
		}
		free(groups);
		free(groupSizes);
	}
#endif

	free(temp);

}

/*
* Find an initial guess for the approximated symmetry direction,
* which in the case of mirror symmetry is vec perpendicular to the mirror plane,
* and in other cases, the symmetry axis.
* This is done using least-mean-squares, which provides 3 guesses, times 3 if we try to remove outliers
*/
void findSymmetryDirection(Molecule *m, double  ***dirs, int *n_dirs, OperationType type) {
	int* groupSizes = (int*)malloc(sizeof(int) * m->groupNum());
	double** groupAverages = (double**)malloc(sizeof(double*) * m->groupNum());
	int *outliers = (int*)malloc(sizeof(int) * m->groupNum());
	int i, j;
	double **testDir;
	double median;
	double zero[3] = { 0.0, 0.0, 0.0 };
	std::vector<double *> results;
	int useOrthogonal = true;
	const double A = 2.0; // What is this? This was previously in the options but was always set to 2.

	testDir = (double**)malloc(sizeof(double*) * 3);
	for (i = 0; i < 3; i++) {
		testDir[i] = (double*)malloc(sizeof(double) * 3);
	}

	for (i = 0; i < m->groupNum(); i++) {
		groupSizes[i] = 0;
		groupAverages[i] = (double*)malloc(sizeof(double) * 3);
		groupAverages[i][0] = groupAverages[i][1] = groupAverages[i][2] = 0.0;
	}

	for (i = 0; i < m->size(); i++) {
		groupSizes[m->similar(i) - 1]++;
		groupAverages[m->similar(i) - 1][0] += m->pos()[i][0];
		groupAverages[m->similar(i) - 1][1] += m->pos()[i][1];
		groupAverages[m->similar(i) - 1][2] += m->pos()[i][2];
	}

	for (i = 0; i < m->groupNum(); i++) {
		outliers[i] = false;
		groupAverages[i][0] /= groupSizes[i];
		groupAverages[i][1] /= groupSizes[i];
		groupAverages[i][2] /= groupSizes[i];
	}

	if (type == CS) {
		// For CS 			
		// Assuming that all orbit center-of-masses are 
		// found on the plane of reflection - find it			
		// for some reason - we get a few options
		planeFit(groupAverages, m->groupNum(), testDir, outliers);
	}
	else {
		// For CN and SN with N > 2			
		// Assuming that all orbit center-of-masses are found on the axis of symmetry - find it			
		// for some reason - we get a few options
		lineFit(groupAverages, m->groupNum(), testDir, outliers);
	}

	// if there are not enough groups for reliable outlier detection - don't invoke it.
	if (options.detectOutliers && m->groupNum() >= MIN_GROUPS_FOR_OUTLIERS) {
		for (j = 0; j < 3; j++) {
			double **tempDir;

			tempDir = (double**)malloc(sizeof(double*) * 3);
			for (i = 0; i < 3; i++) {
				tempDir[i] = (double*)malloc(sizeof(double) * 3);
			}

			// 1. Find the distance of each point from the line / plane
			// 2. Find the median m
			// 3. for each distance di, if (di / m > A || m / di > A - remove as outlier)
			// 4. recompute line / plane
			double* dists = (double *)malloc(m->groupNum() * sizeof(double));
			for (i = 0; i < m->groupNum(); i++) {
				if (type == CS) {
					dists[i] = fabs(testDir[j][0] * groupAverages[i][0] +
						testDir[j][1] * groupAverages[i][1] +
						testDir[j][2] * groupAverages[i][2]);
				}
				else {
					dists[i] = computeDistanceFromLine(groupAverages[i], zero, testDir[j]);
				}

			}
			median = findMedian(dists, m->groupNum());
			for (i = 0; i < m->groupNum(); i++) {
				if (dists[i] / median > A || dists[i] / median > A) { // TODO: What is this? Why is this checked twice?
					outliers[i] = true;
				}
			}
			if (type == CS) {
				planeFit(groupAverages, m->groupNum(), tempDir, outliers);
			}
			else {
				lineFit(groupAverages, m->groupNum(), tempDir, outliers);
			}
			for (i = 0; i < 3; i++) {
				results.push_back(tempDir[i]);
			}

			free(dists);
			free(tempDir);
		}
		for (i = 0; i < 3; i++) {
			free(testDir[i]);
		}

	}
	else {
		// just copy...
		for (i = 0; i < 3; i++) {
			results.push_back(testDir[i]);
		}

	}

	if (useOrthogonal) {
		// Go over all dirs and add two orthogonal axes
		std::vector<double*> temp = results;
		for (i = 0; i < temp.size(); i++) {
			double* dir = temp[i];
			double* newDir = (double *)malloc(sizeof(double) * 3);
			double* newDir2 = (double *)malloc(sizeof(double) * 3);
			double norm, scal;
			if (fabs(dir[0]) < MINDOUBLE) {
				newDir[0] = 1.0; newDir[1] = 0.0; newDir[2] = 0.0;
				newDir2[0] = 0.0; newDir2[1] = -dir[2]; newDir2[2] = dir[1];
			}
			else if (fabs(dir[1]) < MINDOUBLE) {
				newDir[0] = -dir[2]; newDir[1] = 0.0; newDir[2] = dir[0];
				newDir2[0] = newDir2[2] = 0.0; newDir[1] = 1.0;
			}
			else {
				newDir[0] = -dir[1]; newDir[1] = dir[0]; newDir[2] = 0.0;
				newDir2[0] = 0.0; newDir2[1] = -dir[2]; newDir2[2] = dir[1];
			}
			// normalize first
			norm = sqrt(newDir[0] * newDir[0] + newDir[1] * newDir[1] + newDir[2] * newDir[2]);
			newDir[0] /= norm; newDir[1] /= norm; newDir[2] /= norm;
			// remove projection of first from second (both already orthogonal to dir
			scal = newDir[0] * newDir2[0] + newDir[1] * newDir2[1] + newDir[2] * newDir2[2];
			newDir2[0] -= scal * newDir[0]; newDir2[1] -= scal * newDir[1]; newDir2[2] -= scal * newDir[2];
			// normalize second
			norm = sqrt(newDir2[0] * newDir2[0] + newDir2[1] * newDir2[1] + newDir2[2] * newDir2[2]);
			newDir2[0] /= norm; newDir2[1] /= norm; newDir2[2] /= norm;

			results.push_back(newDir);
			results.push_back(newDir2);

		}
	}



	// initialize results array
	*n_dirs = results.size();
	(*dirs) = (double**)malloc(sizeof(double*) *(*n_dirs));
	for (i = 0; i < *n_dirs; i++) {
		(*dirs)[i] = results[i];
	}

	for (i = 0; i < m->groupNum(); i++) {
		free(groupAverages[i]);
	}
	free(outliers);
	free(groupSizes);
	free(groupAverages);
	free(testDir);
}

void estimatePerm(Molecule* m, int *perm, double *dir, OperationType type) {
	int isImproper = (type != CN) ? true : false;
	int isZeroAngle = (type == CS) ? true : false;
	int maxGroupSize = m->getMaxGroupSize();
	int *group = (int*)malloc(sizeof(int) * maxGroupSize);
	int *used = (int*)malloc(sizeof(int) * m->size());
	int i, j, k, l;
	double rotaionMatrix[3][3];
	double tmpMatrix[3][3] = { { 0.0, -dir[2], dir[1] }, { dir[2], 0.0, -dir[0] }, { -dir[1], dir[0], 0.0 } };
	double angle;
	double **rotated = (double**)malloc(sizeof(double*) * m->size());
	struct distRecord * distances = (struct distRecord *)malloc(sizeof(struct distRecord) * maxGroupSize * maxGroupSize);
	int factor = (isImproper ? (-1) : 1);
	int tableSize = 0;
	int orbitDone, orbitSize, orbitStart;
	int left;
	angle = isZeroAngle ? 0.0 : (2 * M_PI / options.opOrder);

	LOG(debug) << "estimatePerm called";

	for (j = 0; j < m->size(); j++) {
		rotated[j] = (double *)malloc(sizeof(double) * 3);
		rotated[j][0] = rotated[j][1] = rotated[j][2] = 0;
		perm[j] = -1;
	}

	// Prepare the rotation matrix
	for (j = 0; j < 3; j++) {
		for (k = 0; k < 3; k++) {
			rotaionMatrix[j][k] =
				((j == k) ? cos(angle) : 0) +
				(factor - cos(angle)) * dir[j] * dir[k] +
				sin(angle) * tmpMatrix[j][k];
		}
	}

	// Run the operation on the current point set
	for (j = 0; j < m->size(); j++) {
		for (k = 0; k < 3; k++) {
			for (l = 0; l < 3; l++) {
				rotated[j][k] += rotaionMatrix[k][l] * m->pos()[j][l];
			}
		}
		LOG(debug) << boost::format("%d (%4.2f, %4.2f, %4.2f) -> (%4.2f, %4.2f, %4.2f)") %
			j %
			m->pos()[j][0] % m->pos()[j][1] % m->pos()[j][2] %
			rotated[j][0] % rotated[j][1] % rotated[j][2];
	}

	// run over the groups
	for (i = 0; i < m->groupNum(); i++) {
		// Get the group
		int groupSize = m->getGroup(i + 1, group);
		for (j = 0; j < m->size(); j++) {
			used[j] = 0;
		}

		// compute the distance matrix
		for (j = 0; j < groupSize; j++) {
			for (k = 0; k < groupSize; k++) {
				int index = j * groupSize + k;
				distances[index].row = group[j];
				distances[index].col = group[k];
				distances[index].distance = 0;
				for (l = 0; l < 3; l++) {
					distances[index].distance +=
						(m->pos()[group[j]][l] - rotated[group[k]][l]) *
						(m->pos()[group[j]][l] - rotated[group[k]][l]);
				}
				distances[index].distance = sqrt(distances[index].distance);
			}
		}

		tableSize = groupSize * groupSize;

		// Sort the distances			
		qsort(distances, tableSize, sizeof(struct distRecord), distComp);

		stringstream groupstrm;
		for (j = 0; j < groupSize; j++)
			groupstrm << group[j] << " ";
		LOG(debug) << "Working on group " << groupstrm.str();


		left = groupSize;
		// Go over the sorted group, and set the permutation
		for (j = 0; j < tableSize && left > 0; j++) {
			int enoughForFullOrbit = left >= options.opOrder;
			int row = distances[j].row;
			int col = distances[j].col;

			// If we have used this item already - skip it.
			if (perm[row] != -1)
				continue;
			LOG(debug) << row << " " << col << " " << distances[j].distance;

			// If we do not have enought to full groups, set all remaining items to themselves
			if (left == 1 || (type == CN && !enoughForFullOrbit)) {
				for (k = 0; k < groupSize; k++) {
					if (used[group[k]] == 0) {
						perm[group[k]] = group[k];
						LOG(debug) << "set " << group[k] << "<->" << group[k];
					}
				}
				break;
			}
			if (options.opOrder == 2) {
				// Special treatment - only size 1 and two orbits are allowed 
				// If both elements are not yet set, use the element.			
				if (perm[row] == -1 && perm[col] == -1)  {
					perm[row] = col;
					perm[col] = row;
					LOG(debug) << "set " << row << "<->" << col;
					left -= (row == col) ? 1 : 2;
				}
			}
			else {
				// we now want to complete an orbit. 
				if (perm[row] == -1 && used[col] == 0) {
					perm[row] = col;
					used[col] = 1;
					LOG(debug) << "set " << row << "<->" << col;
					left--;
				}
				else {
					continue;
				}

				// if this is an orbit of size one...
				if (row == col) continue;

				// If there is no more room for full orbit, must be SN and size two orbit
				if (type == SN && !enoughForFullOrbit) {
					perm[col] = row;
					used[row] = 1;
					LOG(debug) << "set " << row << "->" << col;
					left--;
					continue;
				}

				// Run until an orbit is complete
				orbitDone = false;
				orbitStart = row;
				orbitSize = 1;

				while (!orbitDone) {
					if (orbitSize == options.opOrder - 1) {
						LOG(debug) << "Closing orbit";
						row = col;
						col = orbitStart;
						orbitDone = true;
					}
					else {
						// Search for the next orbit element
						for (k = j + 1; k < tableSize; k++) {
							if (distances[k].row == col && used[distances[k].col] == 0 &&
								distances[k].col != distances[k].row) {
								if (orbitStart == distances[k].col) {
									if (type == SN && orbitSize == 1) {
										// we have now closed an orbit of size 2			
										orbitDone = true;
									}
									else {
										continue;
									}
								}
								row = distances[k].row;
								col = distances[k].col;
								orbitSize++;
								break;
							}
						}
					}
					perm[row] = col;
					used[col] = 1;
					LOG(debug) << "set " << row << "->" << col;
					left--;
				}
			}
		}

	}

	// verify that the orbits are correct?

	for (j = 0; j < m->size(); j++) {
		free(rotated[j]);
	}

	free(rotated);
	free(distances);
	free(group);
}

void lineFit(double **points, int nPoints, double **dirs, int *outliers) {
	// taken from http://www.mapleprimes.com/forum/linear-regression-in-3d
	double A[3] = { 0, 0, 0 };
	double matrix[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };

	int realNum = 0;

	double norm;
	int i, j, k;

	// Compute A
	for (i = 0; i < nPoints; i++) {
		if (!outliers[i]) {
			A[0] += points[i][0];
			A[1] += points[i][1];
			A[2] += points[i][2];
			realNum++;
		}
	}
	A[0] /= realNum;
	A[1] /= realNum;
	A[2] /= realNum;

	// Compute matrix
	for (i = 0; i < nPoints; i++) {
		if (!outliers[i]) {
			for (j = 0; j < 3; j++) {
				for (k = 0; k < 3; k++) {
					matrix[j][k] += (points[i][j] - A[j]) * (points[i][k] - A[k]);
				}
			}
		}
	}

	vector<EigenResult> eigens = GetEigens(matrix);

	// We just try the three lines?...	
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			dirs[i][j] = eigens[i].vector[j]; // copyMat[j + 1][i + 1];
		}
	}

	for (i = 0; i < 3; i++) {
		norm = sqrt(dirs[i][0] * dirs[i][0] + dirs[i][1] * dirs[i][1] + dirs[i][2] * dirs[i][2]);
		dirs[i][0] /= norm;
		dirs[i][1] /= norm;
		dirs[i][2] /= norm;

	}
}
void planeFit(double **points, int nPoints, double **dirs, int* outliers) {
	// taken from http://www.mapleprimes.com/forum/linear-regression-in-3d
	double A[3] = { 0, 0, 0 };
	double matrix[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };

	double norm;
	int i, j, k;
	int realNum = 0;

	// Compute A
	for (i = 0; i < nPoints; i++) {
		if (!outliers[i]) {
			realNum++;
			A[0] += points[i][0];
			A[1] += points[i][1];
			A[2] += points[i][2];
		}
	}
	A[0] /= realNum;
	A[1] /= realNum;
	A[2] /= realNum;

	// Compute matrix
	for (i = 0; i < nPoints; i++) {
		if (outliers[i]) continue;
		for (j = 0; j < 3; j++) {
			for (k = 0; k < 3; k++) {
				matrix[j][k] += (points[i][j] - A[j]) * (points[i][k] - A[k]);
			}
		}
	}

	// compute the matrix's eigenvalues and eigenvectors.
	vector<EigenResult> eigens = GetEigens(matrix);

	// We just try the three planes...	
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			dirs[i][j] = eigens[i].vector[j]; // copyMat[j + 1][i + 1];
		}
	}
	for (i = 0; i < 3; i++) {
		norm = sqrt(dirs[i][0] * dirs[i][0] + dirs[i][1] * dirs[i][1] + dirs[i][2] * dirs[i][2]);
		dirs[i][0] /= norm;
		dirs[i][1] /= norm;
		dirs[i][2] /= norm;

	}
}

