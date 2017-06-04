import math

import numpy as np

from csm.calculations.approx.approximators import CythonHungarianApproximator, ClassicApproximator, NewApproximator
from csm.calculations.approx.dirs import find_symmetry_directions
from csm.calculations.basic_calculations import process_results, CSMState
from csm.calculations.constants import MINDOUBLE, MAXDOUBLE
from csm.fast import approximate_perm_classic, approximate_perm_hungarian
from csm.fast import CythonPermuter

from csm.calculations.exact_calculations import csm_operation
from csm.molecule.molecule import Molecule


def approx_calculation(op_type, op_order, molecule, sn_max=8, use_best_dir=False, get_orthogonal=True, detect_outliers=False,use_chains=False, hungarian=True, print_approx=False, dirs=None, *args, **kwargs):
    if not dirs:
        dirs = find_symmetry_directions(molecule, use_best_dir, get_orthogonal, detect_outliers, op_type)

    if hungarian:
        approximator_cls = CythonHungarianApproximator
    else:
        approximator_cls = ClassicApproximator
    #approximator_cls=NewApproximator


    if op_type == 'CH':  # Chirality
        #sn_max = op_order
        # First CS
        approximator = approximator_cls('CS', 2, molecule, dirs, print_approx)
        best_result = approximator.find_best_perm()
        if best_result.csm > MINDOUBLE:
            # Try the SN's
            approximator = approximator_cls('SN', 2, molecule, dirs, print_approx)
            for op_order in range(2, sn_max + 1, 2):
                approximator._op_order=op_order
                result = approximator.find_best_perm()
                if result.csm < best_result.csm:
                    best_result = result._replace(op_type = 'SN', op_order = op_order)
                if best_result.csm < MINDOUBLE:
                    break

    else:
        approximator = approximator_cls(op_type, op_order, molecule, dirs, print_approx)
        best_result = approximator.find_best_perm()

    best_result = process_results(best_result)
    #TODO: is calc_local supposed to be supported in approx?
    #if calc_local:
    #    local_csm = compute_local_csm(molecule, best_result.perm, best_result.dir, best_result.op_type,
    #                                  best_result.op_order)
    #    best_result = best_result._replace(local_csm=local_csm)

    return best_result

