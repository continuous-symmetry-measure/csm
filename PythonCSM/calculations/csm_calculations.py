from datetime import datetime
import itertools

import math

import numpy as np
from calculations.constants import MINDOUBLE, MAXDOUBLE
from CPP_wrapper.fast import calc_ref_plane
from molecule.normalizations import de_normalize_coords, normalize_coords
from CPP_wrapper.fast import CythonPermuter, SinglePermPermuter, TruePermChecker,PQPermChecker, CythonPIP
import logging
from recordclass import recordclass

np.set_printoptions(precision=6)

logger = logging.getLogger("csm")

__author__ = 'Itay, Devora, Yael'

CSMState = recordclass('CSMState', ('molecule',
                                    'op_order',
                                    'op_type',
                                    'csm',
                                    'perm',
                                    'dir',
                                    'd_min',
                                    'symmetric_structure',
                                    'local_csm',
                                    'perm_count',))
CSMState.__new__.__defaults__ = (None,) * len(CSMState._fields)

# When this property is set by an outside caller, it is called every permutation iteration with the current CSMState
# This is useful for writing all permutations to file during the calculation
csm_state_tracer_func = None


def process_results(results):
    """
    Final normalizations and de-normalizations
    :param results: CSM old_calculations results
    :param csm_args: CSM args
    """
    #    results.molecule.set_norm_factor(molecule.norm_factor)
    masses = [atom.mass for atom in results.molecule.atoms]
    normalize_coords(results.symmetric_structure, masses)

    results.molecule.de_normalize()
    results.symmetric_structure = de_normalize_coords(results.symmetric_structure, results.molecule.norm_factor)


def perm_count(op_type, op_order, molecule, keep_structure, *args, **kwargs):
    def _len_group_permuter(group_size, cycle_size, add_cycles_of_two):
        """
        Returns the length of the group permuter.
        :param group_Size: Group size
        :param cycle_size: Cycle size
        :param add_cycles_of_two: When true, cycles of size of two are also legal
        :return: Number of permutations that will be returned by the group_permuter
        This is done by enumerating the different *kinds* of cycle_structs, and calculating the number of permutations
        for each kind.
        A *kind* of cycle struct is: (n1 cycles of size c1, n2 cycles of size c2, ...). Where c1..cn are either 1, cycle_size
        or 2 (if add_cycles_of_two is True), and the sum of ni*ci adds to group_size
        """
        def log_fact(n):
            """
            Returns the log(n!)
            """
            result = 0.0
            for i in range(1,n+1):
                result += math.log(i)
            return result

        def log_nCr(n, r):
            """
            Calculates the log of n choose r
            """
            if r > n:
                raise ValueError("r must be <= n for nCr")
            if r < 0:
                raise ValueError("r must be non-negative")
            if n <= 0:
                raise ValueError("n must be positive")

            return log_fact(n) - log_fact(n-r) - log_fact(r)

        def generate_kinds():
            """
            Generates all cycle_struct kinds with the given group_size and cycle_lengths.
            A kind is a mapping from cycle_size to the number of cycles
            """
            add_two = add_cycles_of_two and cycle_size!=2  # Do we add cycles of two?
            for num_cycle_size in range(0, group_size // cycle_size + 1):
                if add_two:
                    left = group_size - cycle_size * num_cycle_size  # amount left
                    for num_two in range(0, left // 2 + 1):
                        yield {cycle_size: num_cycle_size, 2: num_two}
                else:
                    yield {cycle_size: num_cycle_size}

        def log_structs_of_kind(kind):
            """
            Returns the log of the number of cycle_structs of a given kind.
            :param kind: The kind of cycle_struct
            :return: The number of cycle_structs of kind.
            The calculation is done by choosing the cycles.

            For example, if there's only one cycle of size 2, there are group_size choose 2 different structs of this kind.
            If there are two cycles of size 2, there are group_size choose 2 * (group_size-2) choose 2 / 2
            The reasoning here is - first choose the first cycle, then choose the second cycle. Divide by two because the
            order of the two cycles doesn't matter.

            For 3 cycles of size 4 we get:
            group_size choose 4 * (group_size-4) choose 4 * (group_size-8) choose 4 / 3!
            We choose the first group of 4, then the second, then the third, but divide by 3! because the order of groups
            doesn't matter.

            If the group size is 12 and we allow 2 cycles of size 3 and 3 cycles of size 2, we get
            (12C3 * 9C3 / 2) * (6C2 * 4C2 * 2C2) / 2!
            We first choose the groups of size 3, then the groups of size 2.

            The above expressions can be further simplified, but we don't do it - this verison is more readable and performs
            very well.
            """
            log_result = 0
            size_left = group_size
            for cycle_size, cycle_count in kind.items():
                for i in range(cycle_count):
                    log_result += log_nCr(size_left, cycle_size)
                    size_left -= cycle_size
                log_result -= log_fact(cycle_count)
            return log_result

        def log_perms_in_kind(kind):
            """
            Returns the log of the number of permutations of a cycle_struct of a specific kind
            The number of permutations for a cycle of size n is (n-1)! so we just go and add that (add, because
            we're returning the log
            """
            log_result = 0
            for cycle_size, cycle_count in kind.items():
                log_result += cycle_count * log_fact(cycle_size-1)
            return log_result

        count = 0
        for kind in generate_kinds():
            log_num_structs = log_structs_of_kind(kind)
            log_num_perms = log_perms_in_kind(kind)
            count += math.exp(log_num_structs + log_num_perms)
        return count
    if keep_structure:
        permuter = CythonPermuter(molecule, op_order, op_type, PQPermChecker, perm_class=CythonPIP)
        for perm in permuter.permute():
            if permuter.count%1000000==0:
                print("counted", permuter.count, "permutations thus far...")
        count=permuter.count
    else:
        count=1
        groups=molecule.equivalence_classes
        for group in groups:
            count*=_len_group_permuter(len(group), op_order, op_type=='SN')
    return int(count)

def exact_calculation(op_type, op_order, molecule, keep_structure=False, perm=None, calc_local=False, *args, **kwargs):
    if keep_structure:
        permchecker=PQPermChecker
    else:
        permchecker=TruePermChecker
        print("Total perms:", perm_count(op_type, op_order, molecule, keep_structure))

    permuter_class=CythonPermuter

    if op_type == 'CH':  # Chirality
        sn_max = op_order
        # First CS
        best_result = csm_operation('CS', 2, molecule, permuter_class, permchecker, perm)
        best_result.op_type='CS'
        if best_result.csm > MINDOUBLE:
            # Try the SN's
            for op_order in range(2, sn_max + 1, 2):
                result = csm_operation('SN', op_order, molecule, permuter_class, permchecker, perm)
                if result.csm < best_result.csm:
                    best_result = result
                    best_result.op_type='SN'
                    best_result.op_order=op_order
                if best_result.csm > MINDOUBLE:
                    break

    else:
        best_result = csm_operation(op_type, op_order, molecule, permuter_class, permchecker, perm)

    process_results(best_result)
    if calc_local:
        best_result.local_csm = compute_local_csm(molecule, best_result.perm, best_result.dir, best_result.op_type,
                                                  best_result.op_order)

    return best_result


def csm_operation(op_type, op_order, molecule, permuter_class=CythonPermuter, permchecker=TruePermChecker, perm=None):
    """
    Calculates minimal csm, dMin and directional cosines by applying permutations
    that keep the similar atoms within the group.
    Once it finds the optimal permutation , calls the CreateSymmetricStructure on the optimal permutation
    :param current_calc_data: current old_calculations data object
    :param args: The CSM arguments
    :return: A dictionary with all the results: csm, dMin, perm and direction
    """
    start_time=datetime.now()
    best_csm = CSMState(molecule=molecule, op_type=op_type, op_order=op_order, csm=MAXDOUBLE)
    traced_state = CSMState(molecule=molecule, op_type=op_type, op_order=op_order)

    if perm:
        permuter = SinglePermPermuter(np.array(perm), molecule, op_order, op_type)
        logger.debug("SINGLE PERM")
    else:
        permuter = permuter_class(molecule, op_order, op_type, permchecker)

    for calc_state in permuter.permute():
        if permuter.count%1000000==0:
            print("calculated for", permuter.count, "permutations thus far...\t Time:", datetime.now()-start_time)
        csm, dir = calc_ref_plane(op_order, op_type=='CS', calc_state)
        if csm_state_tracer_func:
            traced_state.csm = csm
            traced_state.perm = calc_state.perm
            traced_state.dir = dir
            csm_state_tracer_func(traced_state)

        if csm < best_csm.csm:
            best_csm.csm = csm
            best_csm.dir = dir
            best_csm.perm = calc_state.perm

    best_csm.perm_count=permuter.count
    if best_csm.csm == MAXDOUBLE:
        # failed to find csm value for any permutation
        raise ValueError("Failed to calculate a csm value for %s" % op_type)
    best_csm.d_min = 1.0 - (best_csm.csm / 100 * op_order / (op_order - 1))
    best_csm.symmetric_structure = create_symmetric_structure(molecule, best_csm.perm, best_csm.dir, best_csm.op_type,
                                                              best_csm.op_order, best_csm.d_min)
    return best_csm

def create_rotation_matrix(iOp, op_type, op_order, dir):
    is_improper = op_type != 'CN'
    is_zero_angle = op_type == 'CS'
    W = np.array([[0.0, -dir[2], dir[1]], [dir[2], 0.0, -dir[0]], [-dir[1], dir[0], 0.0]])
    rot = np.zeros((3, 3))
    angle = 0.0 if is_zero_angle else 2 * np.pi * iOp / op_order
    factor = -1 if is_improper and (iOp % 2) == 1 else 1

    # The rotation matrix is calculated similarly to the Rodrigues rotation matrix. The only
    # difference is that the matrix is also a reflection matrix when factor is -1.
    #
    # This is why we took the old C++ code instead of applying the Rodrigues formula directly.
    for s in range(3):
        for t in range(3):
            ang = np.cos(angle) if s == t else 0
            rot[s][t] = ang + ((factor - np.cos(angle)) * dir[s] * dir[t] + np.sin(angle) * W[s][t])

    return rot


def compute_local_csm(molecule, perm, dir, op_type, op_order):
    size = len(molecule.atoms)
    cur_perm = [i for i in range(size)]
    local_csm = np.zeros(size)
    m_pos = np.asarray([np.asarray(atom.pos) for atom in molecule.atoms])

    for i in range(op_order):
        rot = create_rotation_matrix(i, op_type, op_order, dir)

        # set permutation
        cur_perm = [perm[cur_perm[j]] for j in range(size)]

        # apply rotation to each atoms
        rotated = rot @ m_pos[cur_perm[i]]
        difference = rotated - m_pos[j]
        square = np.square(difference)
        sum = np.sum(square)
        local_csm[i] = sum * (100.0 / (2 * op_order))
    return local_csm


def create_symmetric_structure(molecule, perm, dir, op_type, op_order, d_min):
    logger.debug('create_symmetric_structure called')

    cur_perm = np.arange(len(perm))  # array of ints...
    size = len(perm)
    m_pos = np.asarray([np.asarray(atom.pos) for atom in molecule.atoms])
    symmetric = np.copy(m_pos)
    logger.debug("in atoms:")
    logger.debug(symmetric)

    normalization = d_min / op_order

    ########calculate and apply transform matrix#########
    ###for i<OpOrder
    for i in range(1, op_order):
        # get rotation
        rotation_matrix = create_rotation_matrix(i, op_type, op_order, dir)
        logger.debug("Rotation matrix:\n")
        logger.debug(rotation_matrix)
        # rotated_positions = m_pos @ rotation_matrix

        # set permutation
        cur_perm = [perm[cur_perm[j]] for j in range(size)]

        # add correct permuted rotation to atom in outAtoms
        for j in range(len(symmetric)):
            symmetric[j] += rotation_matrix @ m_pos[cur_perm[j]]
        logger.debug("Out atoms")
        logger.debug(symmetric)

    # apply normalization:
    symmetric *= normalization
    logger.debug("normalized out atoms:")
    logger.debug(symmetric)

    return symmetric
