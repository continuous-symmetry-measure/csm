from calculations.csm_calculations_data import CSMCalculationsData
from permutations import molecule_permuter

__author__ = 'YAEL'

MAXDOUBLE = 100000000.0

from CPP_wrapper import csm


def perform_operation(csm_args, data):
    if 'dir' in csm_args:
        result = csm.FindBestPermUsingDir(data)
    else:
        if csm_args['findPerm']:
            result = csm.FindBestPerm(data)
        else:
            result = csm_operation(data, csm_args['opName'])
    return result


def csm_operation(current_calc_data, op_name):
    result_csm = MAXDOUBLE
    dir = []
    optimal_perm = []
    current_calc_data.dir = [0, 0, 0]

    # calculate csm for each valid permutation & remember minimal
    perm_size = len(current_calc_data.molecule.atoms)
    groups = current_calc_data.molecule.equivalence_classes
    opOrder = current_calc_data.opOrder
    isSN = current_calc_data.operationType=='SN'
    for perm in molecule_permuter(perm_size, groups, opOrder, isSN):
        current_calc_data.perm = perm
        current_calc_data = csm.CalcRefPlane(current_calc_data)
        # check, if it's a minimal csm, update dir and optimal perm
        if current_calc_data.csm < result_csm:
            result_csm = current_calc_data.csm
            dir = current_calc_data.dir[:]
            optimal_perm = perm[:]

    if result_csm == MAXDOUBLE:
        # failed to find csm value for any permutation
        raise ValueError("Failed to calculate a csm value for %s" % op_name)

    current_calc_data.dMin = 1.0 - (result_csm / 100 * current_calc_data.opOrder / (current_calc_data.opOrder - 1))
    current_calc_data.perm = optimal_perm
    current_calc_data.dir = dir
    current_calc_data.csm = result_csm

    return csm.CreateSymmetricStructure(current_calc_data)
