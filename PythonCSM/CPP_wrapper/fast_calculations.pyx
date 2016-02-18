import numpy as np
cimport numpy as np

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

cimport cython
@cython.boundscheck(False)
def cross(np.ndarray[DTYPE_t, ndim=1] a, np.ndarray[DTYPE_t, ndim=1] b):
    assert a.dtype == DTYPE and b.dtype == DTYPE
    cdef DTYPE_t a0 = a[0]
    cdef DTYPE_t a1 = a[1]
    cdef DTYPE_t a2 = a[2]
    cdef DTYPE_t b0 = b[0]
    cdef DTYPE_t b1 = b[1]
    cdef DTYPE_t b2 = b[2]
    cdef DTYPE_t x = a1*b2 - a2*b1
    cdef DTYPE_t y = a2*b0 - a0*b2
    cdef DTYPE_t z = a0*b1 - a1*b0

    return np.array([x,y,z], dtype=DTYPE).T
