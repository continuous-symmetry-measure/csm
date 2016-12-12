from datetime import datetime
import itertools

import math

import numpy as np
from csm.calculations.basic_calculations import compute_local_csm, process_results, CSMState
from csm.calculations.constants import MINDOUBLE, MAXDOUBLE
from csm.fast import calc_ref_plane

from csm.fast import CythonPermuter, SinglePermPermuter
from csm.calculations.permuters import ConstraintPermuter
import logging

np.set_printoptions(precision=6)

logger = logging.getLogger("csm")

__author__ = 'Itay, Devora, Yael'



# When this property is set by an outside caller, it is called every permutation iteration with the current CSMState
# This is useful for writing all permutations to file during the calculation
csm_state_tracer_func = None




def perm_count(op_type, op_order, molecule, keep_structure, print_perms=False, *args, **kwargs):
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


    if not print_perms:
        if not keep_structure:
            count=1
            groups=molecule.equivalence_classes
            for group in groups:
                count*=_len_group_permuter(len(group), op_order, op_type=='SN')
            return int(count)


    traced_state = CSMState(molecule=molecule, op_type=op_type, op_order=op_order)
    permuter = CythonPermuter(molecule, op_order, op_type, keep_structure, precalculate=False)

    for state in permuter.permute():
        if csm_state_tracer_func:
            traced_state=traced_state._replace(csm = '', perm = state.perm, dir = '')
            csm_state_tracer_func(traced_state)
        if permuter.count%1000000==0:
            print("counted", int(permuter.count/1000000), "million permutations thus far...")
        count=permuter.count
    return count



def exact_calculation(op_type, op_order, molecule, sn_max=8, keep_structure=False, perm=None, calc_local=False, no_constraint=False, *args, **kwargs):
    if op_type == 'CH':  # Chirality
        #sn_max = op_order
        # First CS
        best_result = csm_operation('CS', 2, molecule, keep_structure, perm, no_constraint)
        best_result = best_result._replace(op_type='CS') #unclear why this line isn't redundant
        if best_result.csm > MINDOUBLE:
            # Try the SN's
            for op_order in range(2, sn_max + 1, 2):
                result = csm_operation('SN', op_order, molecule, keep_structure, perm, no_constraint)
                if result.csm < best_result.csm:
                    best_result = result._replace(op_type = 'SN', op_order = op_order)
                if best_result.csm < MINDOUBLE:
                    break

    else:
        best_result = csm_operation(op_type, op_order, molecule, keep_structure, perm, no_constraint)

    best_result = process_results(best_result)
    if calc_local:
        local_csm = compute_local_csm(molecule, best_result.perm, best_result.dir, best_result.op_type,
                                      best_result.op_order)
        best_result = best_result._replace(local_csm=local_csm)

    return best_result


def csm_operation(op_type, op_order, molecule, keep_structure=False, perm=None, no_constraint=False):
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
    else:
        permuter = ConstraintPermuter(molecule, op_order, op_type, keep_structure)
        if no_constraint:
            permuter = CythonPermuter(molecule, op_order, op_type, keep_structure)


    for calc_state in permuter.permute():
        if permuter.count%1000000==0:
            print("calculated for", int(permuter.count / 1000000), "million permutations thus far...\t Time:",
                  datetime.now() - start_time)
        csm, dir = calc_ref_plane(op_order, op_type=='CS', calc_state)

        if csm_state_tracer_func:
            traced_state = traced_state._replace(csm=csm, perm=calc_state.perm, dir=dir)
            csm_state_tracer_func(traced_state)

        if csm < best_csm.csm:
            best_csm = best_csm._replace(csm=csm, dir=dir, perm=list(calc_state.perm))



    if best_csm.csm == MAXDOUBLE:
        # failed to find csm value for any permutation
        raise ValueError("Failed to calculate a csm value for %s %d" % (op_type, op_order))

    if not perm:
        print("Number of permutations: %5.4g" % permuter.count)
        print("Number of branches in permutation tree: %5.4g" % permuter.truecount)
        print("Number of dead ends: %5.4g" % permuter.falsecount)


    best_csm = best_csm._replace(perm_count=permuter.count)
    return best_csm


