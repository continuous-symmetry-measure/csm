__author__ = 'zmbq'

cdef extern from "../libcsm/csm.h":
    int RunCSM(int argc, char **argv)
    int SayHello()