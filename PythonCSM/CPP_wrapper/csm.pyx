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


cdef init_csm_data(csmlib.csm_calculation_data &data, csm_args):
    data.molecule = cppize_molecule(csm_args['molecule'])
    data.outAtoms = []
    if 'dir' in csm_args:
        data.dir = csm_args['dir']
    else:
        data.dir = []
    data.csm = 0
    data.dMin = 0
    if 'perm' in csm_args:
        data.perm = csm_args['perm']
    else:
        data.perm = []
    data.localCSM = []
    data.operationType = cs(csm_args['type'])

cdef parse_csm_data(csmlib.csm_calculation_data &data):
    results = {}
    return results


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

def SetCSMOptions(args):
    cdef csmlib.python_cpp_bridge options
    init_options(options, args)

    csmlib.SetCSMOptions(options)

def Calculate():
    output = csmlib.RunCSM()
    result = parse_output(output)
    return result

def TotalNumberOfPemrutations():
    cdef double num
    num = csmlib.TotalNumberOfPermutations()
    return num

def RunSinglePerm(csm_args):
    cdef csmlib.csm_calculation_data data
    init_csm_data(data, csm_args)
    cdef csmlib.csm_calculation_data result = csmlib.RunSinglePerm(data)
    return parse_csm_data(result)
