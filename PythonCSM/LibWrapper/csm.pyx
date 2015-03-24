__author__ = 'zmbq'

#cimport csmlib

from cpython.string cimport PyString_AsString

#cdef char ** to_cstring_array(list_str):
#    cdef char **ret = <char **>malloc(len(list_str) * sizeof(char *))
#    for i in xrange(len(list_str)):
#        ret[i] = PyString_AsString(list_str[i])
#    return ret

#cdef run(args):
#    argv = to_cstring_array(args)
#    result = csmlib.RunCSM(len(args), argv)
#    # free(argv)
#    return result

#cdef hello():
#    return csmlib.SayHello()

cdef hello_python():
    print("Hello from Cython")
    return -17
