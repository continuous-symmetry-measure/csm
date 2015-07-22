""" The Python wrapper of csmlib """

from calculations.molecule import Atom
from calculations.molecule import Molecule
include "misc.pxi"
include "molecule.pxi"

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

    options.displayPerms = args['displayPerms']

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
