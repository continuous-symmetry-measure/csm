import math
from calculations.csm_calculations_data import CSMCalculationsData
from permutations import molecule_permuter

__author__ = 'YAEL'

MAXDOUBLE = 100000000.0

from CPP_wrapper import csm


def perform_operation(csm_args, data):
    """
    Performs the csm calculation according to the options specified in csm_args
    :param csm_args: the csm arguments dictionary
    :param data: current calculations data object
    :return: the calculations data object with the permutation, csm, dMin and direction values updated
    """
    if 'dir' in csm_args:
        result = csm.FindBestPermUsingDir(data)
    else:
        if csm_args['findPerm']:
            result = csm.FindBestPerm(data)
        else:
            result = csm_operation(data, csm_args['opName'])
    return result


def csm_operation(current_calc_data, op_name):
    """
    Calculates minimal csm, dMin and directional cosines by applying permutations
    that keep the similar atoms within the group.
    Once it finds the optimal permutation , calls the CreateSymmetricStructure on the optimal permutation
    :param current_calc_data: current calculations data object
    :param op_name: operation name
    :return: the calculations data object with the permutation, csm, dMin and direction values updated
    """
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


def total_number_of_permutations(csm_args):
    """
    Calculates the number of permutations to be checked in the csm calculation
    :param csm_args: the csm arguments dictionary
    :return: the number of permutations
    """

    def num_permutations(molecule, operation_order, operation_type):
        total = 1.0
        if operation_order > 2 and operation_type == 'SN':
            # In this case - we enumerate over groups of 1, 2, N
            for group in molecule.equivalence_classes:
                group_size = len(group)
                temp = 0
                fact = math.factorial(group_size)
                # Enumerate over the number of groups of size two
                for k in range(int(group_size / 2) + 1):
                    # Enumerate over the number of groups of size operation_order
                    for j in range(int((group_size - 2 * k) / operation_order) + 1):
                        temp += fact / (math.factorial(group_size - j * operation_order - k * 2) *
                                        math.pow(operation_order, j) * math.pow(2, k) * math.factorial(k))
                total *= temp
        else:
            for group in molecule.equivalence_classes:
                group_size = len(group)
                temp = 0
                fact = math.factorial(group_size)
                for j in range(int(group_size / operation_order) + 1):
                    temp += fact / (math.factorial(group_size - j * operation_order) *
                                    math.pow(operation_order, j) * math.factorial(j))
                total *= temp
        return total

    if csm_args['type'] != 'CH':
        return num_permutations(csm_args['molecule'], csm_args['opOrder'], csm_args['type'])
    else:
        num_perms = num_permutations(csm_args['molecule'], 2, 'SN')
        for i in range(2, csm_args['sn_max'] + 1, 2):
            num_perms += num_permutations(csm_args['molecule'], i, 'SN')
        return num_perms


