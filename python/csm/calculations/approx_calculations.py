import numpy as np
import math
import logging
import munkres
from csm.calculations.constants import MINDOUBLE, MAXDOUBLE
from csm.fast import CythonPermuter, SinglePermPermuter
from csm.fast import estimate_perm #as cython_estimate_perm
from csm.fast import external_get_eigens as cppeigen
from csm.calculations.exact_calculations import csm_operation
from csm.calculations.basic_calculations import CSMState, process_results
from csm.molecule.molecule import Molecule


logger = logging.getLogger("csm")

def approx_calculation(op_type, op_order, molecule, detect_outliers=False,use_chains=False, hungarian=False, print_approx=False, dir=None, *args, **kwargs):
    results= find_best_perm(op_type, op_order, molecule, detect_outliers, use_chains, hungarian, print_approx, dir)
    #print(is_legal(results.perm, molecule))
    return process_results(results)

def trivial_calculation(op_type, op_order, molecule, use_chains=True, *args, **kwargs):
    if molecule.chains and use_chains:
        best = CSMState(molecule=molecule, op_type=op_type, op_order=op_order, csm=MAXDOUBLE)
        chainkeys=list(molecule.chains.keys())
        chain_permutations = []
        dummy = Molecule.dummy_molecule(len(molecule._chains), molecule.chain_equivalences)
        permuter = CythonPermuter(dummy, op_order, op_type, keep_structure=False, precalculate=False)
        for state in permuter.permute():
            chain_permutations.append(list(state.perm))
        for chainperm in chain_permutations:
            perm = [-1 for i in range(len(molecule))]
            for f_index in range(len(chainperm)):
                f_chain_key=chainkeys[f_index]
                t_chain_key=chainkeys[chainperm[f_index]]
                f_chain=molecule.chains[f_chain_key]
                t_chain = molecule.chains[t_chain_key]
                for i in range(len(f_chain)):
                    perm[f_chain[i]]=t_chain[i]

            result = csm_operation(op_type, op_order, molecule, keep_structure=False, perm=perm)
            if result.csm < best.csm:
                best = result

    else:
        perm = [i for i in range(len(molecule))]
        best = csm_operation(op_type, op_order, molecule, keep_structure=False, perm=perm)

    return process_results(best)

def is_legal(perm, molecule):
    legal=True

    switch=True
    for key in molecule.chains:
        if 0 in molecule.chains[key] and perm[0] in molecule.chains[key]:
            switch=False
        break


    switched=True
    for i in range(len(perm)):
        if i in molecule.chains[key] and perm[i] in molecule.chains[key]:
            switched=False
        else:
            switched=True
        if switched != switch:
            test=perm[i]
            first= molecule.atoms[i].chain
            second= molecule.atoms[test].chain
            legal=False
            break

    return legal


def find_best_perm(op_type, op_order, molecule, detect_outliers, use_chains, hungarian, print_approx, dir):
    chain_permutations = []
    dummy = Molecule.dummy_molecule(len(molecule.chains), molecule.chain_equivalences)
    permuter = CythonPermuter(dummy, op_order, op_type, keep_structure=False, precalculate=False)
    for state in permuter.permute():
        chain_permutations.append([i for i in state.perm])

    best = CSMState(molecule=molecule, op_type=op_type, op_order=op_order, csm=MAXDOUBLE)

    if dir:
        for chainperm in chain_permutations:
            perm = estimate_perm(op_type, op_order, molecule, dir, chainperm, False, hungarian)
            best_for_chain_perm = csm_operation(op_type, op_order, molecule, keep_structure=False, perm=perm)
            if best_for_chain_perm.csm < best.csm:
                best = best_for_chain_perm

    elif op_type == 'CI' or (op_type == 'SN' and op_order == 2):
        # if inversion:
        # not necessary to calculate dir, use geometrical center of structure
        dir = [1.0, 0.0, 0.0]
        for chainperm in chain_permutations:
            perm = estimate_perm(op_type, op_order, molecule, dir, chainperm, False, hungarian)
            best_for_chain_perm = csm_operation(op_type, op_order, molecule, keep_structure=False, perm=perm)
            if best_for_chain_perm.csm < best.csm:
                best = best_for_chain_perm


    else:
        dirs = find_symmetry_directions(molecule, detect_outliers, op_type)
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


                    if interim_results.csm < best_for_chain_perm.csm:
                        diff = best_for_chain_perm.csm - interim_results.csm
                        best_for_chain_perm = interim_results
                        if best_for_chain_perm.csm < best.csm:
                            best = best_for_chain_perm

    return best


min_group_for_outliers = 10



def find_symmetry_directions(molecule, detect_outliers, op_type):
    # get average position of each equivalence group:
    group_averages = []
    for group in molecule.equivalence_classes:
        sum = [0, 0, 0]
        for index in group:
            sum += molecule.Q[index]
        average = sum / len(group)
        group_averages.append(average)
    group_averages = np.array(group_averages)

    dirs = dir_fit(group_averages)
    if detect_outliers and len(molecule.equivalence_classes) > min_group_for_outliers:
        dirs = dirs_without_outliers(dirs, group_averages, op_type)
    dirs = dirs_orthogonal(dirs)
    return dirs


def normalize_dir(dir):
    dir = np.array(dir)  # lists also get sent to this function, but can't be /='d
    norm = math.sqrt(math.fabs(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]))
    dir /= norm
    return dir

def dir_fit(positions):
    # get vector of the average position
    sum = np.einsum('ij->j', positions)
    average = sum / len(positions)

    mat = np.zeros((3, 3), dtype=np.float64, order="c")
    for pos in positions:
        for j in range(3):
            for k in range(3):
                mat[j][k] += (pos[j] - average[j]) * (pos[k] - average[k])

    # computer eigenvalues and eigenvectors
    # lambdas - list of 3 eigenvalues of matrix (function demands them, but we won't be using them here)
    # dirs- list of 3 eigenvectors of matrix
    dirs = np.zeros((3, 3), dtype=np.float64, order="c")
    lambdas = np.zeros((3), dtype=np.float64, order="c")
    cppeigen(mat, dirs, lambdas)

    # normalize result:
    dirs = np.array([normalize_dir(dir) for dir in dirs])

    return dirs

def dirs_without_outliers(dirs, positions, op_type):
    def compute_distance_from_line(group_avg_point, test_dir_end):
        def magnitude(point1, point2):
            vec = point2 - point1
            return math.sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2])

        def get_U(Point, LineStart, LineEnd, LineMag):
            return (((Point[0] - LineStart[0]) * (LineEnd[0] - LineStart[0])) +
                    ((Point[1] - LineStart[1]) * (LineEnd[1] - LineStart[1])) +
                    ((Point[2] - LineStart[2]) * (LineEnd[2] - LineStart[2]))) / (LineMag * LineMag)

        def get_intersect(LineStart, LineEnd, U):
            Intersection = [0, 0, 0]
            Intersection[0] = LineStart[0] + U * (LineEnd[0] - LineStart[0]);
            Intersection[1] = LineStart[1] + U * (LineEnd[1] - LineStart[1]);
            Intersection[2] = LineStart[2] + U * (LineEnd[2] - LineStart[2]);
            return Intersection

        linestart = [0, 0, 0]
        linemag = magnitude(test_dir_end, linestart)
        U = get_U(group_avg_point, linestart, test_dir_end, linemag)
        intersection = get_intersect(linestart, test_dir_end, U)
        return magnitude(group_avg_point, intersection)
    print("======================detecting outliers============================")
    more_dirs = []
    for dir in dirs:
        # 1.Find the distance of each point from the line/plane
        dists=[]
        for i in range(len(positions)):
            pos = positions[i]
            if op_type == 'CS':
                dist = math.fabs(dir[0] * pos[0] + dir[1] * pos[1] + dir[2] * pos[2])
            else:
                dist = compute_distance_from_line(pos, dir)
            dists.append(dist)

        # 2. Find median of the distances m
        median = np.median(np.array(dists))
        # 3. for each distance di, if di > 2*m, remove it as an outlier
        with_outliers_removed = list()
        for i in range(len(positions)):
            pos = positions[i]
            dist = dists[i]
            if not (dist / median > 2 or dist / median > 2):
                with_outliers_removed.append(pos)
        # 4. recompute dirs
        outliers_dirs = dir_fit(with_outliers_removed)
        more_dirs += list(outliers_dirs)
    return np.array(more_dirs)

def dirs_orthogonal(dirs):
    added_dirs = list(dirs)

    for dir in dirs:
        if math.fabs(dir[0]) < MINDOUBLE:
            dir1 = [1.0, 0.0, 0.0]
            dir2 = [0.0, -dir[2], dir[1]]
        elif math.fabs(dir[1]) < MINDOUBLE:
            dir1 = [-dir[2], 0.0, dir[0]]
            dir2 = [0.0, 1.0, 0.0]
        else:
            dir1 = [-dir[1], dir[0], 0.0]
            dir2 = [0.0, -dir[2], dir[1]]
        # normalize dir1:
        dir1 = normalize_dir(dir1)
        # remove projection of dir1 from dir2
        scal = dir1[0] * dir2[0] + dir1[1] * dir2[1] + dir1[2] * dir2[2]
        dir2[0] -= scal * dir1[0]
        dir2[1] -= scal * dir1[1]
        dir2[2] -= scal * dir1[2]
        # normalize dir2
        dir2 = normalize_dir(dir2)
        added_dirs.append(dir1)
        added_dirs.append(dir2)

    return np.array(added_dirs)


