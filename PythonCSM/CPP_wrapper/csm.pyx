""" The Python wrapper of csmlib """
__author__ = 'zmbq'

cimport csmlib

def SayHello():
    return csmlib.SayHello()

def RunCSM(args):
    args = [__file__] + args
    encoded = [a.encode('UTF8') for a in args]
    return csmlib.RunCSM(encoded)

