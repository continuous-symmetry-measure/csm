import numpy as np
from csm.calculations.constants import MINDOUBLE, MAXDOUBLE
from csm.fast import CythonPermuter, SinglePermPermuter, TruePermChecker, PQPermChecker, CythonPIP
from csm.fast import external_get_eigens as cppeigen
from csm.calculations.csm_calculations import csm_operation, CSMState, create_rotation_matrix


def find_best_perm(op_type, op_order, molecule, detect_outliers):
    if op_type == 'CI' or (op_type == 'SN' and op_order == 2):
        # if inversion:
        # not necessary to calculate dir, use geometrical center of structure
        dir = [1.0, 0.0, 0.0]
        perm = estimate_perm(op_type, op_order, molecule, dir)
        best = csm_operation(op_type, op_order, molecule, SinglePermPermuter, TruePermChecker, perm)

    else:
        best = CSMState(molecule=molecule, op_type=op_type, op_order=op_order, csm=MAXDOUBLE)
        dirs = find_symmetry_directions(molecule, detect_outliers, op_type)
        for dir in dirs:
            # find permutation for this direction of the symmetry axis
            perm = estimate_perm(op_type, op_order, molecule, dir)

            # solve using this perm until it converges:
            interim_results = csm_operation(op_type, op_order, molecule, SinglePermPermuter, TruePermChecker, perm)
            best_for_this_dir = interim_results
            # (can add in conditions that interim_results.csm>1e-4 and interim_results.csm - old_resutls.csm>0.01 later)
            # (can also use max_iters instead of fixed 50, but 50 is whats given in the article and we don't care right now)
            for i in range(50):
                perm = estimate_perm(op_type, op_order, molecule, interim_results.dir)
                interim_results = csm_operation(op_type, op_order, molecule, SinglePermPermuter, TruePermChecker, perm)
                if interim_results.csm < best_for_this_dir.csm:
                    best_for_this_dir = interim_results

            if best_for_this_dir.csm < best.csm:
                best = best_for_this_dir

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
    additional_dirs = list(dirs)
    for dir in dirs:
        dists = []
        # 1.Find the distance of each point from the line/plane
        for pos in positions:
            if op_type == 'CS':
                dists.append(math.abs(dir[0] * pos[0] + dir[1] * pos[1] + dir[2] * pos[2]))
            else:
                dists.append(np.linalg.norm(dir - pos))

        # 2. Find median of the distances m
        median = np.median(np.array(dists))
        # 3. for each distance di, if di > 2*m, remove it as an outlier
        with_outliers_removed = list()
        for i in range(len(positions)):
            if dists[i] < 2 * median:
                with_outliers_removed.append(positions[i])
        # 4. recompute dirs
        additional_dirs = additional_dirs + list(dir_fit(with_outliers_removed))

    return additional_dirs


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


def estimate_perm(op_type, op_order, molecule, dir):
    # step two: first estimation of the symmetry measure: first permutation
    # apply the symmetry operation once, using symm_element from step one

    # create rotation matrix
    rotation_mat = create_rotation_matrix(1, op_type, op_order, dir)

    # run rotation matrix on atoms
    rotated = (rotation_mat @ molecule.Q.T).T

    # create permutation:
    perm = [-1] * len(molecule)

    #np_distances = np.ones((len(molecule), len(molecule))) * MAXDOUBLE


    # permutation creation is done by group:
    for group in molecule.equivalence_classes:
        if len(group) == 1:
            perm[group[0]] = group[0]
            #np_distances[group[0]][group[0]]=0

        else:
            # measure all the distances between all the points in X to all the points in Y
            distances={}
            for i in group:
                for j in group:
                    a = rotated[i]
                    b = molecule.Q[j]
                    distances[(i,j)]= np.linalg.norm(a - b)
                    #np_distances[i][j] = np.linalg.norm(a - b)

            # recursive:
            # find a pair of points (X,Y) minimally distant from each other (eg, minimum of distances_matrix
            # new_matrix = remove x row, y column from distances_matrix

            # INEFFICIENT METHOD:

            while True:
                min_distance = MAXDOUBLE
                for (i, j) in distances:
                    if distances[(i, j)] < min_distance:
                        min_distance = distances[(i, j)]
                        pair = (i, j)
                i = pair[0]
                j = pair[1]
                for k in group:
                    for t in group:
                        try:
                            if k == i:
                                del distances[(k, t)]
                            elif t == j:
                                del distances[(k, t)]
                        except:
                            pass

                perm[i] = j
                if min_distance == MAXDOUBLE:
                    break

    #np_distances permutation is done by matrix:
    #view= np_distances[:,:]



    # return paired sets, as permutation
    return perm
