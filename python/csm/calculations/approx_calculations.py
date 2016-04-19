import numpy as np
import math
from csm.calculations.constants import MINDOUBLE, MAXDOUBLE
from csm.fast import CythonPermuter, SinglePermPermuter, TruePermChecker, PQPermChecker, CythonPIP
from csm.fast import external_get_eigens as cppeigen
from csm.calculations.csm_calculations import csm_operation, CSMState, create_rotation_matrix

def approx_calculation(op_type, op_order, molecule, detect_outliers=True, *args, **kwargs):
    return find_best_perm(op_type, op_order, molecule, detect_outliers)

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
            best_for_this_dir = old_results= interim_results
            i=0
            max_iterations=50
            while (i<max_iterations and math.fabs(old_results.csm- interim_results.csm)>0.01 and interim_results.csm>0.0001):
                old_results = interim_results
                i+=1
                perm = estimate_perm(op_type, op_order, molecule, interim_results.dir)
                interim_results = csm_operation(op_type, op_order, molecule, SinglePermPermuter, TruePermChecker, perm)
                if interim_results.csm < best_for_this_dir.csm:
                    best_for_this_dir = interim_results


            if best_for_this_dir.csm < best.csm:
                best = best_for_this_dir

            print("best csm is", best.csm)

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

    # permutation creation is done by group:
    for group in molecule.equivalence_classes:
        if len(group) == 1:
            perm[group[0]] = group[0]

        else:
            # measure all the distances between all the points in X to all the points in Y
            size=len(group)
            distances=np.zeros((size,size))
            for i in range(len(group)):
                for j in range(len(group)):
                    a = rotated[group[i]]
                    b = molecule.Q[group[j]]
                    distances[i][j]= np.linalg.norm(a - b)

            # recursive:
            # find a pair of points (X,Y) minimally distant from each other (eg, minimum of distances_matrix
            # new_matrix = remove x row, y column from distances_matrix

            # INEFFICIENT METHOD:

            min_distance = MAXDOUBLE
            maxvals = np.ones(len(group)) * MAXDOUBLE
            left=len(group)
            while left>0:
                if left<op_order:
                    for index in left:
                        perm[index]=index
                        left-=1
                else:
                    (i,j) = (int(np.argmin(distances) / len(distances[0])), int(np.argmin(distances) % len(distances[0])))
                    perm[group[j]] = group[i]
                    distances[:, j] = maxvals
                    distances[i, :] = maxvals
                    left-=1

                     #starting cycle
                    orbitstart=i
                    orbitdone=False
                    orbitlength=1
                    while orbitlength<op_order:
                        (i, j) = (int(np.argmin(distances) / len(distances[0])), int(np.argmin(distances) % len(distances[0])))
                        perm[group[j]] = group[i]
                        distances[:, j] = maxvals
                        distances[i, :] = maxvals
                        left -= 1
                        orbitlength+=1
                    perm[group[orbitstart]]=group[j] #close loop
    return perm




'''
            for i in range(len(distances)): #row
                if left <op_order:
                    j=i
                else:
                    j= np.argmin(distances[i]) #find minimum value for this index
                    if op_order==2 or i ==j: #complete 1 or 2 length cycle
                        perm[group[j]] = group[i] #this line is redundant in the i==j case, but I'm assuming for now the efficiency cost is negligble

                    else: #we're trying to build an orbit
                        1
                #set perm
                perm[group[i]] = group[j]
                distances[:, j] = maxvals #"zero" out the column of j so it cant be selected again
'''

    #np_distances permutation is done by matrix:
    #view= np_distances[:,:]



    # return paired sets, as permutation

