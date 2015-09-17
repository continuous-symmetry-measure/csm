from permutations.lengths import len_molecule_permuter
from permutations.permuters import molecule_permuter
# from CPP_wrapper.permutations import molecule_permuter

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
            result = csm_operation(data, csm_args) # csm_args['opName'], csm_args['molecule'].chains_perms)
    return result


def csm_operation(current_calc_data, csm_args): # op_name, chains_perms):
    """
    Calculates minimal csm, dMin and directional cosines by applying permutations
    that keep the similar atoms within the group.
    Once it finds the optimal permutation , calls the CreateSymmetricStructure on the optimal permutation
    :param current_calc_data: current calculations data object
    :param args: The CSM arguments
    :return: the calculations data object with the permutation, csm, dMin and direction values updated
    """
    op_name = csm_args['opName']
    chains_perms = csm_args['molecule'].chains_perms

    result_csm = MAXDOUBLE
    dir = []
    optimal_perm = []
    current_calc_data.dir = [0, 0, 0]
    # calculate csm for each valid permutation & remember minimal
    for perm in molecule_permuter(len(current_calc_data.molecule.atoms), current_calc_data.molecule.equivalence_classes,
                                  current_calc_data.opOrder, current_calc_data.operationType == 'SN'):
        if csm_args['printPermutations']:
            print(perm)
        current_calc_data.perm = perm
        current_calc_data = csm.CalcRefPlane(current_calc_data)
        # check, if it's a minimal csm, update dir and optimal perm
        if current_calc_data.csm < result_csm:
            (result_csm, dir, optimal_perm) = (current_calc_data.csm, current_calc_data.dir[:], perm[:])

        if chains_perms:
            new_perm = [-1 for i in perm]
            for chain_perm in chains_perms:
                # Apply chain_perm on perm
                for i in range(len(perm)):
                    new_perm[chain_perm[i]] = perm[i]

                current_calc_data.perm = new_perm
                current_calc_data = csm.CalcRefPlane(current_calc_data)
                # check, if it's a minimal csm, update dir and optimal perm
                if current_calc_data.csm < result_csm:
                    (result_csm, dir, optimal_perm) = (current_calc_data.csm, current_calc_data.dir[:], new_perm[:])

    if result_csm == MAXDOUBLE:
        # failed to find csm value for any permutation
        raise ValueError("Failed to calculate a csm value for %s" % op_name)

    current_calc_data.dMin = 1.0 - (result_csm / 100 * current_calc_data.opOrder / (current_calc_data.opOrder - 1))
    current_calc_data.perm = optimal_perm
    current_calc_data.dir = dir
    current_calc_data.csm = result_csm

    return csm.CreateSymmetricStructure(current_calc_data)


def total_number_of_permutations(csm_args):
    if csm_args['type'] != 'CH':
        return len_molecule_permuter(csm_args['molecule'], csm_args['opOrder'], csm_args['type'])
    else:
        num_perms = len_molecule_permuter(csm_args['molecule'], 2, 'SN')
        for i in range(2, csm_args['sn_max'] + 1, 2):
            num_perms += len_molecule_permuter(csm_args['molecule'], i, 'SN')
        return num_perms
