import math

import numpy as np
cimport numpy as np
cimport csmlib

from calculations.constants import MAXDOUBLE, ZERO_IM_PART_MAX

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

ITYPE = np.int
ctypedef np.int_t ITYPE_t

cimport cython


cdef _calc_B(np.ndarray[DTYPE_t, ndim=1, mode="c"] B, int size, np.ndarray[DTYPE_t, ndim=2, mode="c"]Q, np.ndarray[DTYPE_t, ndim=2, mode="c"]Q_, double sintheta):
    csmlib.calc_B(<double *>B.data, size, <double (*)[3]>Q.data, <double (*)[3]>Q_.data, sintheta)


@cython.boundscheck(False)
def calc_B(np.ndarray[DTYPE_t, ndim=1, mode="c"] B not None, int size, np.ndarray[DTYPE_t, ndim=2, mode="c"]Q not None, np.ndarray[DTYPE_t, ndim=2, mode="c"]Q_ not None, double sintheta):
    _calc_B(B, size, Q, Q_, sintheta)


def PolynomialRoots(coeffs):
    cdef double coeffs_v[7]
    cdef double zeror[7]
    cdef double zeroi[7]

    cdef int i
    for i in range(7):
        coeffs_v[i] = coeffs[i]

    csmlib.rpoly(coeffs_v, 6, zeror, zeroi)
    result = []
    for i in range(7):
        result.append(complex(zeror[i], zeroi[i]))

    return result

