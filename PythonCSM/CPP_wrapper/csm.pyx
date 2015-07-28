""" The Python wrapper of csmlib """

include "misc.pxi"
include "molecule.pxi"

cdef init_options(csmlib.python_cpp_bridge &options, args):
    options.opType = cs(args['type'])
    options.opName = cs(args['opName'])
    options.opOrder = args['opOrder']

    options.printLocal = args['printLocal']
    options.writeOpenu = args['writeOpenu']

    options.findPerm = args['findPerm']
    options.limitRun = args['limitRun']
    options.timeOnly = args['timeOnly']
    if args['sn_max']:
        options.sn_max = args['sn_max']

    options.detectOutliers = args['detectOutliers']
    options.babelTest = args['babelTest']

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

    options.molecule = cppize_molecule(args['molecule'])

cdef parse_output(csmlib.csm_output &output):
    cdef int i;
    results = {}

    results['molecule'] = pythonize_molecule(output.molecule)

    outAtoms = []
    for i in range(output.molecule.atoms.size()):
        outAtoms.append(vector_double_to_tuple(output.outAtoms[i]))
    results['outAtoms'] = outAtoms

    results['csm'] = output.csm
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
