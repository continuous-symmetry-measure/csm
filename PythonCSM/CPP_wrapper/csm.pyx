""" The Python wrapper of csmlib """

import sys
cimport csmlib

def SayHello():
    return csmlib.SayHello()

def convert_string(s):
    return s.encode('UTF8')

def RunCSM(arg_dict):
    options = csmlib.csm_options

    options.findPerm = arg_dict['findPerm']
    options.format = convert_string(arg_dict['format'])
    options.inFile = arg_dict['inFile'].fileno()

    # Fill the rest of the options
    return csmlib.RunCSM(options)

