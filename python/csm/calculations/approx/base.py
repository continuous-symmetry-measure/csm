"""
A base class for all Approximators - implementations of approximate algorithms
"""
import math
import numpy as np

from csm.calculations.approx.dirs import find_symmetry_directions
from csm.calculations.basic_calculations import CSMState
from csm.calculations.constants import MAXDOUBLE
from csm.calculations.exact_calculations import csm_operation
from csm.molecule.molecule import Molecule
from csm.fast import CythonPermuter


class Approximator:
    def __init__(self, op_type, op_order, molecule, use_best_dir=False, get_orthogonal=True, detect_outliers=False,use_chains=False, print_approx=False, dirs=None):
        self._op_type = op_type
        self._op_order = op_order
        self._molecule = molecule
        self._print_approx = print_approx
        if dirs:
            self._initial_directions=dirs
        else:
            self._initial_directions=self._choose_initial_directions(molecule, use_best_dir, get_orthogonal, detect_outliers,
                                                            op_type)

    def _print(self, *strings):
        if self._print_approx:
            print(*strings)

    def approximate(self):
        # the basic steps of direction-based approximation are as follows:
        # 0. precalculation of any variables that only need to be calculated once
        self._precalculate()
        best = CSMState(molecule=self._molecule, op_type=self._op_type, op_order=self._op_order, csm=MAXDOUBLE)
        # 1. choose an initial direction
        for dir in self._initial_directions:
            #calculate on the basis of that permutation as detailed in the function
            result= self._approximate_from_initial_dir(dir)
            # 5. repeat from 1, using a different starting direction (assuming more than one)
            if result.csm<best.csm:
                best=result
        # 6. return the best result
        return best

    def _precalculate(self):
        pass

    def _choose_initial_directions(self, molecule, use_best_dir, get_orthogonal, detect_outliers,
                                                            op_type):
        # if inversion:
        # not necessary to calculate dir, use geometrical center of structure
        if self._op_type == 'CI' or (self._op_type == 'SN' and self._op_order == 2):
            initial_directions= [[1.0, 0.0, 0.0]]

        else:
            initial_directions = find_symmetry_directions(molecule, use_best_dir, get_orthogonal, detect_outliers,
                                                            op_type)

        self._print("There are", len(initial_directions), "initial directions to search for the best permutation")
        return initial_directions

    def _approximate_from_initial_dir(self, dir):
        # 2. using the direction, create a permutation
        # 3. using the permutation and it's results, choose another direction
        # 4. repeat 2-3 until an endpoint
        raise NotImplementedError
