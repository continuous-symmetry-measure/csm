from libcpp.vector cimport vector
from libcpp.string cimport string

# Wrapper of the cplusplus library

cdef extern from "interface.h":
    void HelloWorld()
    int Add(int a, int b)
    int AddList(const vector[int])
    string Print(const string message)
    string GetName(const string first, const string last)
    string Concat(const vector[string] parts)
