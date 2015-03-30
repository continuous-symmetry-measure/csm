__author__ = 'zmbq'

# Wrapper of the cplusplus library

cdef extern from "interface.h":
    void HelloWorld()
    int Add(int a, int b)
