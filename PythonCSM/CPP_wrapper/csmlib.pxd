from libcpp.vector cimport vector
from libcpp.string cimport string
__author__ = 'zmbq'

cdef extern from "csmlib.h":
    int SayHello();
    int RunCSM(const vector[string] args);

