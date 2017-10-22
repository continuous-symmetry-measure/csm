import math

import numpy as np

from csm.calculations.approx.approximators import HungarianApproximator, OldApproximator, ManyChainsApproximator
from csm.calculations.approx.dirs import DirectionChooser
from csm.calculations.basic_calculations import process_results, CSMState
from csm.calculations.constants import MINDOUBLE, MAXDOUBLE
from csm.calculations.exact_calculations import Calculation, PlaceHolderOperation


class ApproxCalculation(Calculation):
    def __init__(self, operation, molecule, approx_algorithm='hungarian', sn_max=8, use_best_dir=False, get_orthogonal=True, detect_outliers=False, print_approx=False, dirs=None, *args, **kwargs):
        super().__init__(operation, molecule)
        if approx_algorithm == 'hungarian':
            self.approximator_cls = HungarianApproximator
        if approx_algorithm == 'greedy':
            self.approximator_cls = OldApproximator
        if approx_algorithm == 'many-chains':
            self.approximator_cls = ManyChainsApproximator
        self.direction_chooser=DirectionChooser(molecule, operation.type, operation.order, use_best_dir, get_orthogonal, detect_outliers, dirs)
        self.sn_max=sn_max
        self.calc()

    def calc(self):
        op_type=self.operation.type
        op_order=self.operation.order
        molecule=self.molecule

        # step two: run the appropriate approximator
        if op_type == 'CH':  # Chirality
            # sn_max = op_order
            # First CS
            approximator = self.approximator_cls('CS', 2, molecule, self.direction_chooser,  print_approx=False)
            best_result = approximator.approximate()
            if best_result.csm > MINDOUBLE:
                # Try the SN's
                approximator = self.approximator_cls('SN', 2, molecule, self.direction_chooser, print_approx=False)
                for op_order in range(2, self.sn_max + 1, 2):
                    approximator._op_order = op_order
                    result = approximator.approximate()
                    if result.csm < best_result.csm:
                        best_result = result._replace(op_type='SN', op_order=op_order)
                    if best_result.csm < MINDOUBLE:
                        break
        else:
            approximator = self.approximator_cls(op_type, op_order, molecule, self.direction_chooser, print_approx=False)
            best_result = approximator.approximate()
        # step three: process and return results
        self._csm_result= process_results(best_result)
        return self.result

    def log(self, *args, **kwargs):
        pass



def approx_calculation(op_type, op_order, molecule, approx_algorithm='hungarian', sn_max=8, use_best_dir=False, get_orthogonal=True, detect_outliers=False, print_approx=False, dirs=None, *args, **kwargs):
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
    ac=ApproxCalculation(PlaceHolderOperation(op_type, op_order), molecule, approx_algorithm, sn_max, use_best_dir, get_orthogonal, detect_outliers, print_approx, dirs)
    return ac.result



