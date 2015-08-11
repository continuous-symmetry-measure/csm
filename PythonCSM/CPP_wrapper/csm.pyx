""" The Python wrapper of csmlib """

include "misc.pxi"
include "molecule.pxi"

from calculations.csm_calculations_data import CSMCalculationsData

cdef init_options(csmlib.python_cpp_bridge &options, args):
    options.opType = cs(args['type'])
    options.opName = cs(args['opName'])
    options.opOrder = args['opOrder']

    options.writeOpenu = args['writeOpenu']

    if args['sn_max']:
        options.sn_max = args['sn_max']

    options.detectOutliers = args['detectOutliers']

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


cdef python_data_obj_to_csm_data(csmlib.csm_calculation_data &data, python_data_object):
    data.molecule = cppize_molecule(python_data_object.molecule)
    data.outAtoms = python_data_object.outAtoms
    data.dir = python_data_object.dir
    data.csm = python_data_object.csm
    data.dMin = python_data_object.dMin
    data.perm = python_data_object.perm
    data.localCSM = python_data_object.localCSM
    data.operationType = cs(python_data_object.operationType)
    data.chMinOrder = python_data_object.chMinOrder
    data.chMinType = cs(python_data_object.chMinType)


cdef parse_csm_data(csmlib.csm_calculation_data &data):
    result = CSMCalculationsData()
    result.molecule = pythonize_molecule(data.molecule)
    outAtoms = []
    for i in range(data.molecule.atoms.size()):
        outAtoms.append(vector_double_to_tuple(data.outAtoms[i]))
    result.outAtoms = outAtoms

    result.csm = data.csm
    result.dir = vector_double_to_tuple(data.dir)
    result.dMin = data.dMin
    result.localCSM = vector_double_to_list(data.localCSM)
    result.perm = vector_int_to_list(data.perm)
    result.chMinOrder = data.chMinOrder
    result.chMinType = data.chMinType

    return result

def SetCSMOptions(args):
    cdef csmlib.python_cpp_bridge options
    init_options(options, args)

    csmlib.SetCSMOptions(options)

def TotalNumberOfPemrutations():
    cdef double num
    num = csmlib.TotalNumberOfPermutations()
    return num

def RunSinglePerm(python_data_obj):
    cdef csmlib.csm_calculation_data data
    python_data_obj_to_csm_data(data, python_data_obj)
    cdef csmlib.csm_calculation_data result = csmlib.RunSinglePerm(data)
    return parse_csm_data(result)

def FindBestPermUsingDir (python_data_obj):
    cdef csmlib.csm_calculation_data data
    python_data_obj_to_csm_data(data, python_data_obj)
    cdef csmlib.csm_calculation_data result = csmlib.FindBestPermUsingDir(data)
    return parse_csm_data(result)

def FindBestPerm (python_data_obj):
    cdef csmlib.csm_calculation_data data
    python_data_obj_to_csm_data(data, python_data_obj)
    cdef csmlib.csm_calculation_data result = csmlib.FindBestPerm(data)
    return parse_csm_data(result)

def CsmOperation (python_data_obj):
    cdef csmlib.csm_calculation_data data
    python_data_obj_to_csm_data(data, python_data_obj)
    cdef csmlib.csm_calculation_data result = csmlib.CsmOperation(data)
    return parse_csm_data(result)

def ComputeLocalCSM  (python_data_obj):
    cdef csmlib.csm_calculation_data data
    python_data_obj_to_csm_data(data, python_data_obj)
    cdef csmlib.csm_calculation_data result = csmlib.ComputeLocalCSM(data)
    return parse_csm_data(result)

def DisplayPermutations():
    csmlib.DisplayPermutations()

def GetPermutations(size, groupSize, addGroupsOfTwo):
    cdef vector[vector[int]] c_perms
    c_perms = csmlib.GetPermutations(size, groupSize, addGroupsOfTwo)

    perms = []
    cdef i
    for i in range(c_perms.size()):
        perm = []
        for j in range(size):
            perm.append(c_perms[i][j])
        perms.append(perm)

    return perms
