from libcpp.vector cimport vector

def cs(s):
    """ Converts a Python string to a C++ string """
    return s.encode('UTF8')

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
