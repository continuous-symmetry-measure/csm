"""
A base class for all Approximators - implementations of approximate algorithms
"""
import math
import numpy as np

from csm.calculations.basic_calculations import CSMState
from csm.calculations.constants import MAXDOUBLE
from csm.calculations.exact_calculations import csm_operation
from csm.molecule.molecule import Molecule
from csm.fast import CythonPermuter

class Approximator:
    def __init__(self, op_type, op_order, molecule, dirs, print_approx):
        self._op_type = op_type
        self._op_order = op_order
        self._molecule = molecule
        self._print_approx = print_approx
        self._dirs = dirs
        self._chain_permutations = self._calc_chain_permutations()

    def _calc_chain_permutations(self):
        chain_permutations = []
        dummy = Molecule.dummy_molecule_from_size(len(self._molecule.chains), self._molecule.chain_equivalences)
        permuter = CythonPermuter(dummy, self._op_order, self._op_type, keep_structure=False, precalculate=False)
        for state in permuter.permute():
            chain_permutations.append([i for i in state.perm])

        return chain_permutations

    def _for_inversion(self, best):
        # if inversion:
        # not necessary to calculate dir, use geometrical center of structure
        dir = [1.0, 0.0, 0.0]
        if self._print_approx:
            if self._op_type == 'SN':
                op_msg = 'S2'
            else:
                op_msg = 'CI'
            print("Operation %s - using just one direction: %s" % (op_msg, dir))

        for chainperm in self._chain_permutations:
            self._print("Calculating for chain permutation ", chainperm)
            perm = self._approximate(dir, chainperm)
            best_for_chain_perm = csm_operation(self._op_type, self._op_order, self._molecule, keep_structure=False,
                                                perm=perm)
            if best_for_chain_perm.csm < best.csm:
                best = best_for_chain_perm

        return best

    def _print(self, *strings):
        if self._print_approx:
            print(*strings)

    def find_best_perm(self):
        self._calc_chain_permutations()
        best = CSMState(molecule=self._molecule, op_type=self._op_type, op_order=self._op_order, csm=MAXDOUBLE)

        if self._op_type == 'CI' or (self._op_type == 'SN' and self._op_order == 2):
            return self._for_inversion(best)

        #else:
        self._print("There are", len(self._dirs), "initial directions to search for the best permutation")
        for dir in self._dirs:
            self._print("Calculating for initial direction: ", dir)
            for chainperm in self._chain_permutations:
                self._print("\tCalculating for chain permutation ", chainperm)
                # find permutation for this direction of the symmetry axis
                perm = self._approximate(dir, chainperm)
                # solve using this perm until it converges:
                old_results = CSMState(molecule=self._molecule, op_type=self._op_type, op_order=self._op_order, csm=MAXDOUBLE)
                best_for_chain_perm = interim_results = csm_operation(self._op_type, self._op_order, self._molecule,
                                                                      keep_structure=False, perm=perm)
                self._print("\t\tfound initial permutation")
                self._print("\t\tfirst pass yielded dir", interim_results.dir,
                          "and CSM " + str(round(interim_results.csm, 5)))
                    # print(perm)

                if best_for_chain_perm.csm < best.csm:
                    best = best_for_chain_perm

                # iterations:
                i = 0
                max_iterations = 50
                while (i < max_iterations and
                           (math.fabs(old_results.csm - interim_results.csm) / math.fabs(
                               old_results.csm) > 0.01 and interim_results.csm < old_results.csm) and interim_results.csm > 0.0001):
                    old_results = interim_results
                    i += 1
                    perm = self._approximate(interim_results.dir, chainperm)
                    interim_results = csm_operation(self._op_type, self._op_order, self._molecule, keep_structure=False,
                                                    perm=perm)

                    self._print("\t\titeration", i, ":")
                    self._print("\t\t\tfound a permutation using dir", old_results.dir, "...")
                    self._print("\t\t\tthere are",
                              len(perm) - np.sum(np.array(perm) == np.array(old_results.perm)),
                              "differences between new permutation and previous permutation")
                    self._print("\t\t\tusing new permutation, found new direction", interim_results.dir)
                    self._print("\t\t\tthe distance between the new direction and the previous direction is:",
                              str(round(np.linalg.norm(interim_results.dir - old_results.dir), 8)))
                    self._print("\t\t\tthe csm found is:", str(round(interim_results.csm, 8)))
                        # print(perm)

                    if interim_results.csm < best_for_chain_perm.csm:
                        diff = best_for_chain_perm.csm - interim_results.csm
                        best_for_chain_perm = interim_results
                        if best_for_chain_perm.csm < best.csm:
                            best = best_for_chain_perm

        return best
