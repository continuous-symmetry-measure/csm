""" The Python wrapper of csmlib """

from libcpp.vector cimport vector
cimport csmlib

def SayHello():
    return csmlib.SayHello()

def cs(s):
    """ Converts a Python string to a C++ string """
    return s.encode('UTF8')

cdef fill_molecule(csmlib.python_cpp_bridge options, args):
    cdef vector[csmlib.python_atom] atoms
    cdef csmlib.python_atom bridge_atom

    for atom in args['molecule']:
        bridge_atom.symbol = cs(atom.symbol)
        bridge_atom.adjacent = atom.adjacent
        bridge_atom.pos = atom.pos
        bridge_atom.mass = atom.mass
        atoms.push_back(bridge_atom)

    options.molecule = atoms

def RunCSM(args):
    cdef csmlib.python_cpp_bridge options

    options.opType = cs(args['type'])
    options.opName = cs(args['opName'])
    options.opOrder = args['opOrder']

    options.printNorm = args['printNorm']
    options.printLocal = args['printLocal']
    options.writeOpenu = args['writeOpenu']

    if args['format']:
        options.format = cs(args['format'])

    options.ignoreHy = args['ignoreHy']
    options.removeHy = args['removeHy']
    options.ignoreSym = args['ignoreSym']

    options.findPerm = args['findPerm']
    options.useMass = args['useMass']
    options.limitRun = args['limitRun']
    options.babelBond = args['babelBond']
    options.timeOnly = args['timeOnly']
    if args['sn_max']:
        options.sn_max = args['sn_max']

    options.detectOutliers = args['detectOutliers']
    options.babelTest = args['babelTest']
    options.keepCenter = args['keepCenter']

    if args['logFileName']:
        options.logFilename = cs(args['logFileName'])

    options.inFilename = cs(args['inFileName'])
    options.fdIn = args['inFile'].fileno()

    options.outFilename = cs(args['outFileName'])
    options.fdOut = args['outFile'].fileno()

    if 'perm' in args:
        options.perm = args['perm']
    else:
        options.perm = []

    if 'dir' in args:
        options.dir = args['dir']
    else:
        options.dir = []

    fill_molecule(options, args)

    print("Calling C++ from Python")
    result = csmlib.RunCSM(options)
    print("Returning to Python")
    return result
