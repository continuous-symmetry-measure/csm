from libcpp.vector cimport vector
from cpython cimport array

import numpy as np
cimport numpy as np

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

def cs(s):
    """ Converts a Python string to a C++ string """
    return s.encode('UTF8')

def ps(cs):
    """ Converts a C++ string to a Python string """
    return cs.decode('UTF8')

cdef vector_double_to_list(const vector[double] &vec):
    cdef int i;

    lst = []
    for i in range(vec.size()):
        lst.append(vec[i])
    return lst

cdef vector_int_to_list(const vector[int] &vec):
    cdef int i;

    lst = []
    for i in range(vec.size()):
        lst.append(vec[i])
    return lst

cdef vector_double_to_tuple(const vector[double] &vec):
    return tuple(vector_double_to_list(vec))

cdef vector[int] list_to_vector_int(lst):
    cdef vector[int] vec;

    for x in lst:
        vec.push_back(x)
    return vec

cdef array.array _int_array_template = array.array('i', [])
cdef inline int *ptr(array.array int_array):
    return int_array.data.as_ints

cdef array.array int_array(int size, zeros=False):
    return array.clone(_int_array_template, size, zeros)
