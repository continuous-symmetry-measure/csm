import datetime
import itertools

import math

import numpy as np
from csm.calculations.basic_calculations import process_results, CSMState
from csm.calculations.constants import MINDOUBLE, MAXDOUBLE, start_time
from csm.fast import calc_ref_plane

from csm.fast import CythonPermuter, SinglePermPermuter
from csm.calculations.permuters import ConstraintPermuter
import logging

from csm.input_output.formatters import format_perm_count

np.set_printoptions(precision=6)


__author__ = 'Itay, Devora, Yael'



# When this property is set by an outside caller, it is called every permutation iteration with the current CSMState
# This is useful for writing all permutations to file during the calculation
csm_state_tracer_func = None

class CSMValueError(ValueError):
    def __init__(self, arg1, CSMState):
        self.arg1 = arg1
        self.CSMState = CSMState
        super().__init__(arg1)


def exact_calculation(op_type, op_order, molecule, sn_max=8, keep_structure=False, perm=None, calc_local=False, no_constraint=False, suppress_print=False, timeout=300, *args, **kwargs):
    """
    Calculates the exact CSM of a molecule by measure the symmetry of every possible permutation of its atoms and choosing the best
    :param op_type: type of symmetry (CS, CN, CH, CI, SN)
    :param op_order: order of symmetry (2, 3, 4...)
    :param molecule: instance of Molecule class whose symmetry is being measured
    :param sn_max: for chirality, the maximum SN symmetry to measure
    :param keep_structure: when True, only permutations that do not break the bond structure of the molecule are measured
    :param perm: default None. If a perm (list of indices) is provided, returns the exact symmetry measure of only that permutation
    :param calc_local: 
    :param no_constraint: when True, the old Cython 
    :param suppress_print: 
    :param timeout: 
    :return: 
    """
    if op_type == 'CH':  # Chirality
        #sn_max = op_order
        # First CS
        best_result = csm_operation('CS', 2, molecule, keep_structure, perm, no_constraint, suppress_print, timeout)
        best_result = best_result._replace(op_type='CS') #unclear why this line isn't redundant
        if best_result.csm > MINDOUBLE:
            # Try the SN's
            for op_order in range(2, sn_max + 1, 2):
                result = csm_operation('SN', op_order, molecule, keep_structure, perm, no_constraint, suppress_print, timeout)
                if result.csm < best_result.csm:
                    best_result = result._replace(op_type = 'SN', op_order = op_order)
                if best_result.csm < MINDOUBLE:
                    break

    else:
        best_result = csm_operation(op_type, op_order, molecule, keep_structure, perm, no_constraint, suppress_print, timeout)

    best_result = process_results(best_result, calc_local)
    return best_result


def csm_operation(op_type, op_order, molecule, keep_structure=False, perm=None, no_constraint=False, suppress_print=False, timeout=300):
    """
    Calculates minimal csm, dMin and directional cosines by applying permutations
    that keep the similar atoms within the group.
    Once it finds the optimal permutation , calls the CreateSymmetricStructure on the optimal permutation
    :param current_calc_data: current old_calculations data object
    :param args: The CSM arguments
    :return: A dictionary with all the results: csm, dMin, perm and direction
    """
    best_csm = CSMState(molecule=molecule, op_type=op_type, op_order=op_order, csm=MAXDOUBLE)
    traced_state = CSMState(molecule=molecule, op_type=op_type, op_order=op_order)

    if perm:
        permuter = SinglePermPermuter(np.array(perm), molecule, op_order, op_type)
    else:
        permuter = ConstraintPermuter(molecule, op_order, op_type, keep_structure, timeout=timeout)
        if no_constraint:
            permuter = CythonPermuter(molecule, op_order, op_type, keep_structure, timeout=timeout)


    for calc_state in permuter.permute():
        if permuter.count%1000000==0:
            print("calculated for", int(permuter.count / 1000000), "million permutations thus far...\t Time:",
                   datetime.datetime.now()-start_time)
        csm, dir = calc_ref_plane(op_order, op_type=='CS', calc_state)

        if csm_state_tracer_func:
            traced_state = traced_state._replace(csm=csm, perm=calc_state.perm, dir=dir)
            csm_state_tracer_func(traced_state)

        if csm < best_csm.csm:
            best_csm = best_csm._replace(csm=csm, dir=dir, perm=list(calc_state.perm))

    if best_csm.csm == MAXDOUBLE:
        # failed to find csm value for any permutation
        best_csm = best_csm._replace(csm=csm, dir=dir, perm=list(calc_state.perm))
        raise CSMValueError("Failed to calculate a csm value for %s %d" % (op_type, op_order), best_csm)

    if not perm and not suppress_print:
        print("Number of permutations: %s" % format_perm_count(permuter.count))
        print("Number of branches in permutation tree: %s" % format_perm_count(permuter.truecount))
        print("Number of dead ends: %s" % format_perm_count(permuter.falsecount))


    best_csm = best_csm._replace(perm_count=permuter.count)
    return best_csm


