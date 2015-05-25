""" The Python wrapper of csmlib """

from libcpp.vector cimport vector
cimport csmlib
from input_output.molecule import Atom

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

cdef fill_molecule(csmlib.python_cpp_bridge &options, args):
    cdef vector[csmlib.python_atom] atoms
    cdef csmlib.python_atom bridge_atom

    for atom in args['molecule']:
        bridge_atom.symbol = cs(atom.symbol)
        bridge_atom.adjacent = atom.adjacent
        bridge_atom.pos = atom.pos
        bridge_atom.mass = atom.mass
        atoms.push_back(bridge_atom)

    options.molecule = atoms

cdef read_molecule(csmlib.csm_output &output):
    cdef int i;

    atoms = []
    for i in range(output.molecule.size()):
        atom = Atom(output.molecule[i].symbol, vector_double_to_tuple(output.molecule[i].pos), useMass=False)
        atom.adjacent = vector_int_to_list(output.molecule[i].adjacent)
        atom._mass = output.molecule[i].mass

        atoms.append(atom)
    return atoms

cdef init_options(csmlib.python_cpp_bridge &options, args):
    options.opType = cs(args['type'])
    options.opName = cs(args['opName'])
    options.opOrder = args['opOrder']

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

    if 'perm' in args:
        options.perm = args['perm']
    else:
        options.perm = []

    if 'dir' in args:
        options.dir = args['dir']
    else:
        options.dir = []

    fill_molecule(options, args)

cdef parse_output(csmlib.csm_output &output):
    cdef int i;
    results = {}
    results['atoms'] = read_molecule(output)
    results['norm'] = output.norm;
    results['numGroups'] = output.numGroups;

    outAtoms = []
    for i in range(output.molecule.size()):
        outAtoms.append(vector_double_to_tuple(output.outAtoms[i]))
    results['outAtoms'] = outAtoms

    results['csm'] = output.csm;
    results['dir'] = vector_double_to_tuple(output.dir)
    results['dMin'] = output.dMin
    results['localCSM'] = vector_double_to_list(output.localCSM)
    results['chMinOrder'] = output.chMinOrder
    results['chMinType'] = output.chMinType
    results['perm'] = vector_int_to_list(output.perm)

    return results

def RunCSM(args):
    cdef csmlib.python_cpp_bridge options
    init_options(options, args)

    output = csmlib.RunCSM(options)
    result = parse_output(output)
    return result
