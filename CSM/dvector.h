/*
 * A bounded double vector class.
 *
 * dvector is a vector of doubles with a lower and upper range. For example dvector(1,3) is a vector with 3 cells - 1, 2, and 3.
 *
 * dvector is here as a replacement to NRUtil's dvector() function
 *
 * Since a lot of CSM code works with 1-based vectors, recreating NRUtil's horrible dvector() was more sensible
 * than using std::vector and updating all the indices everywhere - undoubtedly introducing many many bugs.
 *
 * Written by Itay Zandbank
 */

#ifndef DVECTOR_H
#define DVECTOR_H

#include <vector>

namespace csm_utils
{

	class dvector
	{
	private:
		std::vector<double> _vec;
		int _lowerBound, _upperBound;
		double *_adjustedBuf;  // _adjustedBuf[index] translates to the right cell in the vector

		void adjustBuf()
		{
			_adjustedBuf = _vec.data() - _lowerBound;
		}

	public:
		dvector(int lowerBound, int upperBound)
		{
			_lowerBound = lowerBound;
			_upperBound = upperBound;
			_vec = std::vector<double>(upperBound - lowerBound + 1);
			adjustBuf();
		}

		dvector(const dvector &other)
		{
			// Explicit copy constructor because adjustedBuf can't be copied as is
			_lowerBound = other._lowerBound;
			_upperBound = other._upperBound;
			_vec = other._vec;
			adjustBuf();
		}

		double &operator[](int index)
		{
			// Use _adjustedBuf and not _vec[index-_lowerBound] because the second is a lot slower
			// on some compilers, and CSM has no known bugs in accessing dvectors.
			return _adjustedBuf[index];
		}
	};

}
#endif