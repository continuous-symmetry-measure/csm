"""
This is where the outside wrapper for the approximate calculation is located and called from
"""

import math

import numpy as np
import sys

from csm.calculations.approx.approximators import HungarianApproximator, GreedyApproximator, ManyChainsApproximator, \
    StructuredApproximator
from csm.calculations.approx.dirs import DirectionChooser
from csm.calculations.basic_calculations import process_results, Operation, Calculation
from csm.calculations.constants import MINDOUBLE, MAXDOUBLE
from csm.input_output.readers import check_perm_validity


class ApproxCalculation(Calculation):
    """
    Runs an approximate algorithm to estimate the csm value, using directions to create permutations iteratively
    """
    def __init__(self, operation, molecule, approx_algorithm='hungarian', use_best_dir=False, get_orthogonal=True, detect_outliers=False, dirs=None, keep_structure=False, timeout=100, fibonacci=False, *args, **kwargs):
        """
        Initializes the ApproxCalculation
        :param operation: instance of Operation class or named tuple, with fields for name and order, that describes the symmetry
        :param molecule: instance of Molecule class on which the described symmetry calculation will be performed
        :param approx_algorithm: a string, can be 'hungarian', 'greedy', or 'many-chains'. Chooses the approximator alogrithm to use
        :param use_best_dir: use only the best direction 
        :param get_orthogonal: get orthogonal direction vectors from the main directions
        :param detect_outliers: detect outliers and use the imrpoved direction vectors
        :param dirs: a list of directions to use as initial dire (will override all other dir-choice options)
        """
        super().__init__(operation, molecule)
        if approx_algorithm == 'hungarian':
            self.approximator_cls = HungarianApproximator
        if approx_algorithm == 'greedy':
            self.approximator_cls = GreedyApproximator
        if approx_algorithm == 'many-chains':
            self.approximator_cls = ManyChainsApproximator
        if keep_structure:
            self.approximator_cls = StructuredApproximator

        self.timeout=timeout
        self.direction_chooser=DirectionChooser(molecule, operation.type, operation.order, use_best_dir, get_orthogonal, detect_outliers, dirs, fibonacci)

    def calculate(self):
        op_type=self.operation.type
        op_order=self.operation.order
        molecule=self.molecule

        # step two: run the appropriate approximator
        if op_type == 'CH':  # Chirality
            # First CS
            approximator = self.approximator_cls('CS', 2, molecule, self.direction_chooser, self.log)
            best_result = approximator.approximate()
            if best_result.csm > MINDOUBLE:
                # Try the SN's
                approximator = self.approximator_cls('SN', 2, molecule, self.direction_chooser, self.log)
                for op_order in range(2, self.operation.order + 1, 2):
                    approximator._op_order = op_order
                    result = approximator.approximate()
                    if result.csm < best_result.csm:
                        best_result = result._replace(op_type='SN', op_order=op_order)
                    if best_result.csm < MINDOUBLE:
                        break
        else:
            approximator = self.approximator_cls(op_type, op_order, molecule, self.direction_chooser, self.log, self.timeout)
            best_result = approximator.approximate()
        # step three: process and return results
        self._csm_result= best_result
        return self.result

    def log(self, *args, **kwargs):
        pass


class PrintApprox(ApproxCalculation):
    def log(self, *args, **kwargs):
        print(*args)

def approx_calculation(op_type, op_order, molecule, approx_algorithm='hungarian', sn_max=8, use_best_dir=False, get_orthogonal=True, detect_outliers=False, print_approx=False, dirs=None, keep_structure=False, *args, **kwargs):
    """
    Runs an approximate algorithm to estimate the csm value, using directions to create permutations iteratively

    :param op_type: type of symmetry (CS, CN, CH, CI, SN)
    :param op_order: order of symmetry (2, 3, 4...)
    :param molecule: instance of Molecule class whose symmetry is being measured
    :param approx_algorithm: string, 'hungarian', 'greedy', or 'many-chains'. many chains is hungarian optimized for many chains.
    :param sn_max: for chirality, the maximum SN symmetry to measure
    :param use_best_dir: use only the best direction 
    :param get_orthogonal: get orthogonal direction vectors from the main directions
    :param detect_outliers: detect outliers and use the imrpoved direction vectors
    :param print_approx: 
    :param dirs: a list of directions to use as initial dire
    :return: CSMResult of approximate calculation
    """
    if print_approx:
        ac=PrintApprox(Operation.placeholder(op_type, op_order, sn_max), molecule, approx_algorithm, use_best_dir, get_orthogonal, detect_outliers, dirs, keep_structure)
    else:
        ac=ApproxCalculation(Operation.placeholder(op_type, op_order, sn_max), molecule, approx_algorithm, use_best_dir, get_orthogonal, detect_outliers, dirs, keep_structure)
    ac.calculate()
    return ac.result



