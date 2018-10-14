import datetime
import itertools

import math

import numpy as np
import sys

from csm.calculations.basic_calculations import check_perm_cycles, now, run_time
from csm.calculations.data_classes import CSMState, Operation, CSMResult
from csm.calculations.constants import MINDOUBLE, MAXDOUBLE
from csm.fast import calc_ref_plane
from csm.fast import CythonPermuter, SinglePermPermuter
import logging

from csm.input_output.formatters import format_perm_count
from csm.input_output.formatters import csm_log as print
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


class ExactStatistics:
    def __init__(self, permuter):
        self._perm_count=permuter.count
        self._truecount=permuter.truecount
        self._falsecount=permuter.falsecount

    def write(self, f=sys.stderr):
        f.write("Number of permutations: %s" % format_perm_count(self.perm_count))
        f.write("Number of branches in permutation tree: %s" % format_perm_count(self.num_branches))
        f.write("Number of dead ends: %s" % format_perm_count(self.dead_ends))

    def to_dict(self):
        return {
            "perm count":self.perm_count,
            "number branches":self.num_branches,
            "dead ends":self.dead_ends
        }

    @property
    def dead_ends(self):
        return self._falsecount

    @property
    def perm_count(self):
        return self._perm_count

    @property
    def num_branches(self):
        return self._truecount

class ExactCalculation:
    def __init__(self, operation, molecule, *args, **kwargs):
        """
        A class for running the exact CSM Algorithm
        :param operation: instance of Operation class or named tuple, with fields for name and order, that describes the symmetry
        :param molecule: instance of Molecule class on which the described symmetry calculation will be performed
        :param keep_structure: boolean, when True only permutations that maintain bond integrity will be measured 
        :param perm: a list of atom indiced describing one permutation of the molecule. Default None- When provided,
         only the provided permutation is measured
        :param no_constraint: boolean, default False, when False the constraints algorithm is used for the permuter, 
        when True the old permuter is used
        :param timeout: default 300, the number of seconds the function will run before timing out
        :param callback_func: default None, this function is called for every single permutation calculated with an argument of a single 
        CSMState, can be used for printing in-progress reports, outputting to an excel, etc.
        """
        self.operation=operation
        self.molecule=molecule

    def calculate(self, timeout=300, *args, **kwargs):
        self.start_time = now()
        op_type=self.operation.type
        op_order=self.operation.order
        molecule=self.molecule


        if op_type == 'CH':  # Chirality
            # sn_max = op_order
            # First CS
            best_result = self.csm_operation('CS', 2, molecule, timeout=timeout)
            best_result = best_result._replace(op_type='CS')  # unclear why this line isn't redundant
            if best_result.csm > MINDOUBLE:
                # Try the SN's
                for op_order in range(2, self.operation.order + 1, 2):
                    result = self.csm_operation('SN', op_order, molecule, timeout=timeout)
                    if result.csm < best_result.csm:
                        best_result = result._replace(op_type='SN', op_order=op_order)
                    if best_result.csm < MINDOUBLE:
                        break

        else:
            best_result = self.csm_operation(op_type, op_order, molecule, timeout=timeout)

        overall_stats = self.statistics.to_dict()
        overall_stats["runtime"]=run_time(self.start_time)
        self._csm_result = CSMResult(best_result, self.operation, overall_stats=overall_stats)
        return self.result

    def csm_operation(self, op_type, op_order, molecule, perm=None, timeout=300):
        """
        Calculates minimal csm, directional cosines by applying permutations that keep the similar atoms within the group.
        :param op_type: cannot be CH.
        :param op_order:
        :param molecule:
        :param keep_structure:
        :param perm:
        :param no_constraint:
        :param suppress_print:
        :param timeout:
        :return:
        """
        best_csm = CSMState(molecule=molecule, op_type=op_type, op_order=op_order, csm=MAXDOUBLE)
        traced_state = CSMState(molecule=molecule, op_type=op_type, op_order=op_order)

        if perm:
            permuter = SinglePermPermuter(np.array(perm, dtype="long"), molecule, op_order, op_type)
        else:
            permuter = CythonPermuter(molecule, op_order, op_type, timeout=timeout)

        for calc_state in permuter.permute():
            if permuter.count % 1000000 == 0:
                print("calculated for", int(permuter.count / 1000000), "million permutations thus far...\t Time:",
                      run_time(self.start_time))
            csm, dir = calc_ref_plane(op_order, op_type == 'CS', calc_state)

            if csm < best_csm.csm:
                best_csm = best_csm._replace(csm=csm, dir=dir, perm=list(calc_state.perm))

        self.statistics=ExactStatistics(permuter)

        if best_csm.csm == MAXDOUBLE:
            # failed to find csm value for any permutation
            best_csm = best_csm._replace(csm=csm, dir=dir, perm=list(calc_state.perm))
            raise CSMValueError("Failed to calculate a csm value for %s %d" % (op_type, op_order), best_csm)
        return best_csm

    @staticmethod
    def exact_calculation_for_approx(operation, molecule, perm):
        ec = ExactCalculation(operation, molecule)
        if operation.type == 'CH':  # Chirality
            best_result = ec.csm_operation('CS', 2, molecule,perm=perm)
            if best_result.csm > MINDOUBLE:
                # Try the SN's
                for op_order in range(2, operation.order + 1, 2):
                    result = ec.csm_operation('SN', op_order, molecule,perm=perm)
                    if result.csm < best_result.csm:
                        best_result = result._replace(op_type='SN', op_order=op_order)
                    if best_result.csm < MINDOUBLE:
                        break
        else:
            best_result=ec.csm_operation(operation.type, operation.order, molecule, perm=perm)

        falsecount, num_invalid, cycle_counts, bad_indices = check_perm_cycles(perm, operation)
        best_result=best_result._replace(num_invalid=num_invalid)
        return best_result

    @property
    def result(self):
        return self._csm_result




