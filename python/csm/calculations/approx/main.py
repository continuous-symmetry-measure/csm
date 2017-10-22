import math

import numpy as np

from csm.calculations.approx.approximators import HungarianApproximator, OldApproximator, ManyChainsApproximator
from csm.calculations.basic_calculations import process_results, CSMState
from csm.calculations.constants import MINDOUBLE, MAXDOUBLE



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

    #step one: choose and create the appropriate Approximator
    if approx_algorithm== 'hungarian':
        approximator_cls = HungarianApproximator
    if approx_algorithm== 'greedy':
        approximator_cls = OldApproximator
    if approx_algorithm== 'many-chains':
        approximator_cls=ManyChainsApproximator

    #step two: run the appropriate approximator
    if op_type == 'CH':  # Chirality
        # sn_max = op_order
        # First CS
        approximator = approximator_cls('CS', 2, molecule, use_best_dir, get_orthogonal, detect_outliers, print_approx, dirs)
        best_result = approximator.approximate()
        if best_result.csm > MINDOUBLE:
            # Try the SN's
            approximator = approximator_cls('SN', 2, molecule, use_best_dir, get_orthogonal, detect_outliers,print_approx, dirs)
            for op_order in range(2, sn_max + 1, 2):
                approximator._op_order = op_order
                result = approximator.approximate()
                if result.csm < best_result.csm:
                    best_result = result._replace(op_type='SN', op_order=op_order)
                if best_result.csm < MINDOUBLE:
                    break
    else:
        approximator = approximator_cls(op_type, op_order, molecule, use_best_dir, get_orthogonal, detect_outliers, print_approx, dirs)
        best_result=approximator.approximate()
    #step three: process and return results
    return process_results(best_result)



