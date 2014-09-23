/*
 * A C++ replacement of NRUtil's dmatrix
 *
 * A double matrix with lower and upper bounds on both its indices
 */

#ifndef DMATRIX_H
#define DMATRIX_H

#include <vector>
#include "dvector.h"

namespace csm
{
	class dmatrix
	{
	private:
		int _rowLowerBound, _rowUpperBound, _colLowerBound, _colUpperBound;

		std::vector<dvector> _rows;
		dvector *_adjustedRows;

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

			_adjustedRows = _rows.data() - _rowLowerBound;
		}

		dvector& operator[](int index)
		{
			//return _rows[index - _rowLowerBound];
			return _adjustedRows[index];
		}
	};
}

#endif