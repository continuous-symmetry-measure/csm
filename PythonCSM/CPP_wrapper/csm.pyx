""" The Python wrapper of csmlib """

import sys
cimport csmlib

def SayHello():
    return csmlib.SayHello()

def RunCSM(args):
    try:
        filename = sys.argv[0]
    except KeyError:
        filename = __file__
    args = [filename] + args
    encoded = [a.encode('UTF8') for a in args]
    return csmlib.RunCSM(encoded)

