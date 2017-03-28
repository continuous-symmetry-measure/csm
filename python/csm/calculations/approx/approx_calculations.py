import math

import numpy as np
from csm.fast import CythonPermuter
from csm.fast import estimate_perm  # as cython_estimate_perm

from csm.calculations.basic_calculations import CSMState
from csm.calculations.constants import MINDOUBLE, MAXDOUBLE
from csm.calculations.exact_calculations import csm_operation
from csm.molecule.molecule import Molecule

def find_best_perm(op_type, op_order, molecule, use_chains, hungarian, print_approx, dirs):
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
            perm = estimate_perm(op_type, op_order, molecule, dir, chainperm, False, hungarian)
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
                perm = estimate_perm(op_type, op_order, molecule, dir, chainperm, use_chains, hungarian)
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
                    perm = estimate_perm(op_type, op_order, molecule, interim_results.dir, chainperm, use_chains, hungarian)
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

