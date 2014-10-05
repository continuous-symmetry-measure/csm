/*
 * A C++ replacement of NRUtil's dmatrix
 *
 * A double matrix with lower and upper bounds on both its indices
 *
 * This is a replacement for NRUtil's horrible dmatrix function. Since CSM uses some 1-based matrices,
 * it was easier to reproduce the behavior, rather than fix all the places that counted on such bizarre indexing.
 *
 * Written by Itay Zandbank
 */

#ifndef DMATRIX_H
#define DMATRIX_H

#include <vector>
#include "dvector.h"

namespace csm_utils
{
	class dmatrix
	{
	private:
		int _rowLowerBound, _rowUpperBound, _colLowerBound, _colUpperBound;

		std::vector<dvector> _rows;
		dvector *_adjustedRows;

		void adjustRows()
		{
			_adjustedRows = _rows.data() - _rowLowerBound;
		}

	public:
		dmatrix(int rowLowerBound, int rowUpperBound, int colLowerBound, int colUpperBound)
		{
			_rowLowerBound = rowLowerBound;
			_rowUpperBound = rowUpperBound;
			_colLowerBound = colLowerBound;
			_colUpperBound = colUpperBound;

			_rows = std::vector<dvector>();
			for (int i = _rowLowerBound; i <= _rowUpperBound; i++)
				_rows.push_back(dvector(_colLowerBound, _colUpperBound));
		}

		dmatrix(const dmatrix& other)
		{
			// Explicit copy constructor because _adjustedRows can't be copied as is.
			_rowLowerBound = other._rowLowerBound;
			_rowUpperBound = other._rowUpperBound;
			_colLowerBound = other._colLowerBound;
			_colUpperBound = other._colUpperBound;
			_rows = other._rows;

			adjustRows();
		}

		dvector& operator[](int index)
		{
			// Use _adjustedBuf and not _rows[index-_lowerBound] because the second is a lot slower
			// on some compilers, and CSM has no known bugs in accessing dvectors.
			return _adjustedRows[index];
		}
	};
}

#endif