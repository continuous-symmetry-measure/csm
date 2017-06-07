import math

import numpy as np

from csm.calculations.approx.approximators import HungarianApproximator, OldApproximator, NewChainsApproximator
from csm.calculations.approx.dirs import find_symmetry_directions
from csm.calculations.basic_calculations import process_results, CSMState
from csm.calculations.constants import MINDOUBLE, MAXDOUBLE
from csm.fast import approximate_perm_classic, approximate_perm_hungarian
from csm.fast import CythonPermuter

from csm.calculations.exact_calculations import csm_operation
from csm.molecule.molecule import Molecule


def approx_calculation(op_type, op_order, molecule, sn_max=8, use_best_dir=False, get_orthogonal=True, detect_outliers=False,use_chains=False, hungarian=True, print_approx=False, dirs=None, *args, **kwargs):

    #step one: choose and create the appropriate Approximator
    if hungarian:
        approximator_cls = HungarianApproximator

    #step two: run the appropriate approximator
    if op_type == 'CH':  # Chirality
        # sn_max = op_order
        # First CS
        approximator = approximator_cls('CS', 2, molecule, use_best_dir, get_orthogonal, detect_outliers,use_chains, print_approx, dirs)
        best_result = approximator.approximate()
        if best_result.csm > MINDOUBLE:
            # Try the SN's
            approximator = approximator_cls('SN', 2, molecule, use_best_dir, get_orthogonal, detect_outliers,use_chains, print_approx, dirs)
            for op_order in range(2, sn_max + 1, 2):
                approximator._op_order = op_order
                result = approximator.find_best_perm()
                if result.csm < best_result.csm:
                    best_result = result._replace(op_type='SN', op_order=op_order)
                if best_result.csm < MINDOUBLE:
                    break
    else:
        approximator = approximator_cls(op_type, op_order, molecule, use_best_dir, get_orthogonal, detect_outliers,use_chains, print_approx, dirs)
        best_result=approximator.approximate()
    #step three: process and return results
    return process_results(best_result)



