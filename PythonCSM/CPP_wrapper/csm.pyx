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

cdef vector[int] list_to_vector_int(lst):
    cdef vector[int] vec;

    for x in lst:
        vec.push_back(x)
    return vec

cdef csmlib.python_molecule cppize_molecule(atoms, equivalence_classes):
    """ Convert the Python structures to a CPP-compatible molecule """
    cdef csmlib.python_molecule molecule
    cdef csmlib.python_atom bridge_atom

    for atom in atoms:
        bridge_atom.symbol = cs(atom.symbol)
        bridge_atom.adjacent = atom.adjacent
        bridge_atom.pos = atom.pos
        bridge_atom.mass = atom.mass
        molecule.atoms.push_back(bridge_atom)

    for ec in equivalence_classes:
        molecule.equivalenceClasses.push_back(list_to_vector_int(ec))

    return molecule

cdef pythonize_molecule(const csmlib.python_molecule &molecule):
    cdef int i;

    atoms = []
    for i in range(molecule.atoms.size()):
        atom = Atom(molecule.atoms[i].symbol, vector_double_to_tuple(molecule.atoms[i].pos), useMass=False)
        atom.adjacent = vector_int_to_list(molecule.atoms[i].adjacent)
        atom._mass = molecule.atoms[i].mass

        atoms.append(atom)

    equivalenceClasses = []
    for i in range(molecule.equivalenceClasses.size()):
        oneClass = vector_int_to_list(molecule.equivalenceClasses[i])
        equivalenceClasses.append(oneClass)
    return (atoms, equivalenceClasses)

cdef init_options(csmlib.python_cpp_bridge &options, args):
    options.opType = cs(args['type'])
    options.opName = cs(args['opName'])
    options.opOrder = args['opOrder']

    options.printLocal = args['printLocal']
    options.writeOpenu = args['writeOpenu']

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

    options.molecule = cppize_molecule(args['molecule'], args['equivalence_classes'])

cdef parse_output(csmlib.csm_output &output):
    cdef int i;
    results = {}

    results['atoms'], results['equivalency'] = pythonize_molecule(output.molecule)
    results['norm'] = output.norm;
    results['numGroups'] = output.numGroups;

    outAtoms = []
    for i in range(output.molecule.atoms.size()):
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

def CallInitSimilarity(atoms):
    cdef csmlib.python_molecule molecule;

    molecule = cppize_molecule(atoms, [])
    csmlib.FillEquivalencyClasses(molecule)
    _, equivalency_classes = pythonize_molecule(molecule)

    return equivalency_classes
