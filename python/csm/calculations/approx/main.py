import math

import numpy as np

from csm.calculations.approx.classic import CythonHungarianApproximator, ClassicApproximator
from csm.calculations.approx.dirs import find_symmetry_directions
from csm.calculations.basic_calculations import process_results, CSMState
from csm.calculations.constants import MINDOUBLE, MAXDOUBLE
from csm.fast import approximate_perm_classic, approximate_perm_hungarian
from csm.fast import CythonPermuter

from csm.calculations.exact_calculations import csm_operation
from csm.molecule.molecule import Molecule


def approx_calculation(op_type, op_order, molecule, sn_max=8, use_best_dir=False, get_orthogonal=True, detect_outliers=False,use_chains=False, hungarian=False, print_approx=False, dirs=None, *args, **kwargs):
    if not dirs:
        dirs = find_symmetry_directions(molecule, use_best_dir, get_orthogonal, detect_outliers, op_type)

    if hungarian:
        approximator_cls = ClassicApproximator
    else:
        approximator_cls = CythonHungarianApproximator

    approximator = approximator_cls(op_type, op_order, molecule, dirs, print_approx)

    if op_type == 'CH':  # Chirality
        #sn_max = op_order
        # First CS
        best_result = _find_best_perm('CS', 2, molecule, print_approx, dirs, approximator)
        best_result = best_result._replace(op_type='CS') #unclear why this line isn't redundant
        if best_result.csm > MINDOUBLE:
            # Try the SN's
            for op_order in range(2, sn_max + 1, 2):
                result = approximator.find_best_perm()
                if result.csm < best_result.csm:
                    best_result = result._replace(op_type = 'SN', op_order = op_order)
                if best_result.csm < MINDOUBLE:
                    break

    else:
        best_result = approximator.find_best_perm()

    best_result = process_results(best_result)
    #TODO: is calc_local supposed to be supported in approx?
    #if calc_local:
    #    local_csm = compute_local_csm(molecule, best_result.perm, best_result.dir, best_result.op_type,
    #                                  best_result.op_order)
    #    best_result = best_result._replace(local_csm=local_csm)

    return best_result


def _find_best_perm(op_type, op_order, molecule, print_approx, dirs, approximator):
    chain_permutations = []
    dummy = Molecule.dummy_molecule_from_size(len(molecule.chains), molecule.chain_equivalences)
    permuter = CythonPermuter(dummy, op_order, op_type, keep_structure=False, precalculate=False)
    for state in permuter.permute():
        chain_permutations.append([i for i in state.perm])

    best = CSMState(molecule=molecule, op_type=op_type, op_order=op_order, csm=MAXDOUBLE)

    if op_type == 'CI' or (op_type == 'SN' and op_order == 2):
        # if inversion:
        # not necessary to calculate dir, use geometrical center of structure
        dir = [1.0, 0.0, 0.0]
        for chainperm in chain_permutations:
            perm = approximator(op_type, op_order, molecule, dir, chainperm)
            best_for_chain_perm = csm_operation(op_type, op_order, molecule, keep_structure=False, perm=perm)
            if best_for_chain_perm.csm < best.csm:
                best = best_for_chain_perm


    else:
        if print_approx:
            print("There are", len(dirs), "initial directions to search for the best permutation")
        for chainperm in chain_permutations:
            if print_approx:
                print("Calculating for chain permutation ", chainperm)
            for dir in dirs:
                if print_approx:
                    print("\tCalculating for initial direction: ", dir)
                # find permutation for this direction of the symmetry axis
                perm = approximator(op_type, op_order, molecule, dir, chainperm)
                # solve using this perm until it converges:
                old_results = CSMState(molecule=molecule, op_type=op_type, op_order=op_order, csm=MAXDOUBLE)
                best_for_chain_perm = interim_results = csm_operation(op_type, op_order, molecule,
                                                                   keep_structure=False, perm=perm)
                if print_approx:
                    print("\t\tfound initial permutation")
                    print("\t\tfirst pass yielded dir", interim_results.dir, "and CSM " + str(round(interim_results.csm, 5)))
                    #print(perm)

                if best_for_chain_perm.csm < best.csm:
                    best = best_for_chain_perm

                #iterations:
                i = 0
                max_iterations = 50
                while (i < max_iterations and
                           (math.fabs(old_results.csm - interim_results.csm) / math.fabs(old_results.csm) > 0.01 and interim_results.csm<old_results.csm) and interim_results.csm > 0.0001):
                    old_results = interim_results
                    i += 1
                    perm = approximator(op_type, op_order, molecule, interim_results.dir, chainperm)
                    interim_results = csm_operation(op_type, op_order, molecule, keep_structure=False, perm=perm)

                    if print_approx:
                        print("\t\titeration", i, ":")
                        print("\t\t\tfound a permutation using dir", old_results.dir, "...")
                        print("\t\t\tthere are", len(perm)-np.sum(np.array(perm)==np.array(old_results.perm)), "differences between new permutation and previous permutation")
                        print("\t\t\tusing new permutation, found new direction", interim_results.dir)
                        print("\t\t\tthe distance between the new direction and the previous direction is:", str(round(np.linalg.norm(interim_results.dir - old_results.dir),8)))
                        print("\t\t\tthe csm found is:", str(round(interim_results.csm, 8)))
                        #print(perm)


                    if interim_results.csm < best_for_chain_perm.csm:
                        diff = best_for_chain_perm.csm - interim_results.csm
                        best_for_chain_perm = interim_results
                        if best_for_chain_perm.csm < best.csm:
                            best = best_for_chain_perm

    return best