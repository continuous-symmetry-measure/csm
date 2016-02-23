import math

import numpy as np
cimport numpy as np
cimport csmlib


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


cdef _get_eigen(np.ndarray[DTYPE_t, ndim=2, mode="c"] A, np.ndarray[DTYPE_t, ndim=2, mode="c"]m, np.ndarray[DTYPE_t, ndim=1, mode="c"]lambdas):
    csmlib.GetEigens(<double (*)[3]>A.data, <double (*)[3]>m.data, <double *>lambdas.data)

def get_eigen(np.ndarray[DTYPE_t, ndim=2, mode="c"] A not None, np.ndarray[DTYPE_t, ndim=2, mode="c"]m not None, np.ndarray[DTYPE_t, ndim=1, mode="c"]lambdas not None):
    _get_eigen(A, m, lambdas)


def calc_A_B(op_order, multiplier, sintheta,np.ndarray[ITYPE_t, ndim=2]perms, size, Q):
    # A is calculated according to formula (17) in the paper
    # B is calculated according to formula (12) in the paper

    A = np.zeros((3, 3,))
    B = np.zeros((3,), dtype=np.float64, order="c")  # Row vector for now

    # compute matrices according to current perm and its powers (the identity does not contribute anyway)
    cdef int i
    for i in range(1, op_order):
        # The i'th power of the permutation
        cur_perm = perms[i]
        # Q_ is Q after applying the i'th permutation on atoms (Q' in the article)
        Q_ = np.array([  Q[p] for p in cur_perm  ],  dtype=np.float64, order="c")  # Q'
        # A_intermediate is calculated according to the formula (5) in the paper, as follows:
        # the cross product of Qi, Q_i plus the cross product of Q_i, Qi, summed.
        A+=multiplier[i] * (Q.T.dot(Q_)+Q_.T.dot(Q))
        #B is the sum of the cross products of Q[k] and Q_[k], times sintheta. it is calculated in c++
        calc_B(B, size, Q, Q_, sintheta[i])

    return A, B.T  # Return B as a column vector