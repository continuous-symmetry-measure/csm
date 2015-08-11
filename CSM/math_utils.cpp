/*
* Some utility mathematical functions.
*
* Extracted from mainRot.cpp during the 2014 reorganization.
* Created by Itay Zandbank
*/

#include "math_utils.h"
#include "options.h"
#include <cmath>
#include "logging.h"

using namespace std;


double Magnitude(double *Point1, double *Point2)
{
	double Vector[3];

	Vector[0] = Point2[0] - Point1[0];
	Vector[1] = Point2[1] - Point1[1];
	Vector[2] = Point2[2] - Point1[2];

	return sqrt(Vector[0] * Vector[0] + Vector[1] * Vector[1] + Vector[2] * Vector[2]);
}

double computeDistanceFromLine(double *Point, double *LineStart, double *LineEnd)
{
	double LineMag;
	double U;
	double Intersection[3];

	LineMag = Magnitude(LineEnd, LineStart);

	U = (((Point[0] - LineStart[0]) * (LineEnd[0] - LineStart[0])) +
		((Point[1] - LineStart[1]) * (LineEnd[1] - LineStart[1])) +
		((Point[2] - LineStart[2]) * (LineEnd[2] - LineStart[2]))) /
		(LineMag * LineMag);

	Intersection[0] = LineStart[0] + U * (LineEnd[0] - LineStart[0]);
	Intersection[1] = LineStart[1] + U * (LineEnd[1] - LineStart[1]);
	Intersection[2] = LineStart[2] + U * (LineEnd[2] - LineStart[2]);

	return Magnitude(Point, Intersection);
}


/*
*  This Quickselect routine is based on the algorithm described in
*  "Numerical recipes in C", Second Edition,
*  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
*  This code by Nicolas Devillard - 1998. Public domain.
*/


#define ELEM_SWAP(a,b) { register double t=(a);(a)=(b);(b)=t; }

double findMedian(double arr[], int n)
{
	int low, high;
	int median;
	int middle, ll, hh;

	low = 0; high = n - 1; median = (low + high) / 2;
	for (;;) {
		if (high <= low) /* One element only */
			return arr[median];

		if (high == low + 1) {  /* Two elements only */
			if (arr[low] > arr[high])
				ELEM_SWAP(arr[low], arr[high]);
			return arr[median];
		}

		/* Find median of low, middle and high items; swap into position low */
		middle = (low + high) / 2;
		if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]);
		if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]);
		if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]);

		/* Swap low item (now in position middle) into position (low+1) */
		ELEM_SWAP(arr[middle], arr[low + 1]);

		/* Nibble from each end towards middle, swapM_PIng items when stuck */
		ll = low + 1;
		hh = high;
		for (;;) {
			do ll++; while (arr[low] > arr[ll]);
			do hh--; while (arr[hh]  > arr[low]);

			if (hh < ll)
				break;

			ELEM_SWAP(arr[ll], arr[hh]);
		}

		/* Swap middle item (in position low) back into correct position */
		ELEM_SWAP(arr[low], arr[hh]);

		/* Re-set active partition */
		if (hh <= median)
			low = ll;
		if (hh >= median)
			high = hh - 1;
	}
}

