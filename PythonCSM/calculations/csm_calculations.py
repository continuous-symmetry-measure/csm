import math
from calculations.csm_calculations_data import CSMCalculationsData
from permutations import molecule_permuter
from permutations.permuters import _get_cycle_structs

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
    for perm in molecule_permuter(len(current_calc_data.molecule.atoms), current_calc_data.molecule.equivalence_classes,
                                  current_calc_data.opOrder, current_calc_data.operationType == 'SN'):
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


def _len_all_circle_permutations(size):
    return math.factorial(size-1)


def _len_all_perms_from_cycle_struct(cycle_struct):
    result = 1
    for cycle in cycle_struct:
        result *= _len_all_circle_permutations(len(cycle))

    return result


def _len_group_permuter(group_size, cycle_size, add_cycles_of_two):
    cycle_lengths = {1, cycle_size}
    if add_cycles_of_two:
        cycle_lengths.add(2)
    result = 0
    for cycle_struct in _get_cycle_structs(group_size, cycle_lengths):
        result += _len_all_perms_from_cycle_struct(cycle_struct)

    return result


def len_molecule_permuter(molecule, op_order, op_type):
    result = 1
    for group in molecule.equivalence_classes:
        result *= _len_group_permuter(len(group), op_order, op_type == 'SN')
    return result


def total_number_of_permutations(csm_args):
    if csm_args['type'] != 'CH':
        return len_molecule_permuter(csm_args['molecule'], csm_args['opOrder'], csm_args['type'])
    else:
        num_perms = len_molecule_permuter(csm_args['molecule'], 2, 'SN')
        for i in range(2, csm_args['sn_max'] + 1, 2):
            num_perms += len_molecule_permuter(csm_args['molecule'], i, 'SN')
        return num_perms


def total_number_of_permutations_CPP(csm_args):
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
