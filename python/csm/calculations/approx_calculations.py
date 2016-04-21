import numpy as np
import math
import logging
from csm.calculations.constants import MINDOUBLE, MAXDOUBLE
from csm.fast import CythonPermuter, SinglePermPermuter, TruePermChecker, PQPermChecker, CythonPIP
from csm.fast import external_get_eigens as cppeigen
from csm.calculations.csm_calculations import csm_operation, CSMState, create_rotation_matrix

logger = logging.getLogger("csm")


def approx_calculation(op_type, op_order, molecule, detect_outliers=False, *args, **kwargs):
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
            best_for_this_dir = interim_results = csm_operation(op_type, op_order, molecule, SinglePermPermuter, TruePermChecker, perm)
            logger.debug("permutation is:" + str(perm))
            old_results = CSMState(molecule=molecule, op_type=op_type, op_order=op_order, csm=MAXDOUBLE)
            i = 0
            max_iterations = 50
            while (i < max_iterations and math.fabs(
                        old_results.csm - interim_results.csm) > 0.01 and interim_results.csm > 0.0001):
                old_results = interim_results
                i += 1
                perm = estimate_perm(op_type, op_order, molecule, interim_results.dir)
                interim_results = csm_operation(op_type, op_order, molecule, SinglePermPermuter, TruePermChecker, perm)
                logger.debug("permutation is:"+str(perm))
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
    logger.debug("averages:"+str(group_averages))
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
    more_dirs=[]
    for dir in dirs:
        logger.debug("original dir was:" + str(dir))
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
                logger.debug("not an outlier:"+str(dists[i]))
                with_outliers_removed.append(positions[i])
        # 4. recompute dirs
        outliers_dirs= dir_fit(with_outliers_removed)
        test=list(outliers_dirs)
        more_dirs+=test
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


class distance_record:
    def __init__(self, row, col, dist):
        self.row = row
        self.distance = dist
        self.col = col
    def __str__(self):
        return "distance" + str(self.row) +"," + str(self.col)+": "+str(self.distance)


def approx_create_rotation_matrix(op_type, op_order, dir):
    is_improper = op_type != 'CN'
    is_zero_angle = op_type == 'CS'
    W = np.array([[0.0, -dir[2], dir[1]], [dir[2], 0.0, -dir[0]], [-dir[1], dir[0], 0.0]])
    logger.debug("W:"+str(W))
    rot = np.zeros((3, 3))
    angle = 0.0 if is_zero_angle else 2 * np.pi / op_order
    factor = -1 if is_improper else 1
    # The rotation matrix is calculated similarly to the Rodrigues rotation matrix. The only
    # difference is that the matrix is also a reflection matrix when factor is -1.
    #
    # This is why we took the old C++ code instead of applying the Rodrigues formula directly.
    for s in range(3):
        for t in range(3):
            ang = np.cos(angle) if s == t else 0
            rot[s][t] = ang + ((factor - np.cos(angle)) * dir[s] * dir[t] + np.sin(angle) * W[s][t])

    return rot


def estimate_perm(op_type, op_order, molecule, dir):
    # this is a naive implementation, which works well with the naive c++ build_perm
    # step two: first estimation of the symmetry measure: first permutation
    # apply the symmetry operation once, using symm_element from step one

    # create rotation matrix
    logger.debug("estimatePerm called on dir: "+str(dir))
    rotation_mat = approx_create_rotation_matrix(op_type, op_order, dir)
    logger.debug("rotation matrix:\n"+ str(rotation_mat))
    # run rotation matrix on atoms
    rotated = (rotation_mat @ molecule.Q.T).T

    # create permutation:
    perm = [-1] * len(molecule)

    # np_distances = np.ones((len(molecule), len(molecule))) * MAXDOUBLE


    # permutation creation is done by group:
    for group in molecule.equivalence_classes:
        if len(group) == 1:
            perm[group[0]] = group[0]
            # np_distances[group[0]][group[0]]=0

        else:
            # measure all the distances between all the points in X to all the points in Y
            distances = list()
            for i in group:
                for j in group:
                    a = rotated[i]
                    b = molecule.Q[j]
                    distance = np.linalg.norm(a - b)
                    distances.append(distance_record(j, i, distance))
                    # np_distances[i][j] = np.linalg.norm(a - b)

            distances = newlist = sorted(distances, key=lambda x: x.distance)
            perm = perm_builder(op_type, op_order, group, distances, perm)

    return perm


def perm_builder(op_type, op_order, group, distances, perm):
    # this is intended to serve as a copy-paste of the (buggy) c++ algorithm
    logger.debug("group:"+ str(group))
    left = len(group)
    used = np.zeros(len(perm))  # array to keep track of which values have already been placed
    # Go over the sorted group, and set the permutation
    for record in distances:
        row = record.row
        col = record.col
        # If we have set or used these indexes already - skip this pair.
        if perm[row] != -1 or used[col] != 0:
            continue
        # If we do not have enought to fill cycles, set all remaining items to themselves
        if left <= 1 or (op_type == 'CN' and left < op_order):
            if perm[row] == -1:
                perm[row] = row
                used[row] = 1
                logger.debug("set [" + str(row) +  "] =" + str (row) + " (distance=" + str(record.distance) + ")")
                left -= 1

        elif op_order == 2:
            # Special treatment - only size  1 and two  orbits are allowed
            # If both elements  are not yet set, use the element.
            if perm[row] == -1 and perm[col] == -1:# and used[row] == used[col] == 0:
                perm[row] = col
                perm[col] = row
                used[row] = used[col] = 1
                logger.debug("set [" + str(row) + "] =" + str(col) + " (distance=" + str(record.distance) + ")")
                logger.debug("set [" + str(col) + "] =" + str(row) + " (distance=" + str(record.distance) + ")")
                if row == col:
                    left -= 1
                else:
                    left -= 2
        else:
            # we now want to complete an orbit
            if perm[row] == -1 and used[col] == 0:
                perm[row] = col
                logger.debug("set [" + str(row) + "] =" + str(col) + " (distance=" + str(record.distance) + ")")
                used[col] = 1
                left -= 1
            else:
                continue
            # if this is an orbit of size one, "cycle" is complete, so...
            if row == col:
                continue
            # If there is no more room for full orbit, must be SN and size two orbit
            if type == 'SN' and left < op_order:
                perm[col] = row
                logger.debug("set [" + str(col) + "] =" + str(row) + " (distance=" + str(record.distance) + ")")
                used[row] = 1
                left -= 1
                continue
                # run until orbit is complete:
            orbit_done = False
            orbit_start = row
            orbit_size = 1
            while not orbit_done:
                if orbit_size == op_order - 1:
                    row = col
                    col = orbit_start
                    orbit_done = True
                else:
                    for next_item in distances:
                        flag=False
                        next_row = next_item.row
                        next_col = next_item.col
                        if next_row == col and used[next_col] == 0 \
                                and next_row != next_col:
                            if next_col == orbit_start:
                                if len(group)>3:
                                    flag=True
                                if type=='SN' and orbit_size==1:
                                        orbit_done=True
                                else:
                                        continue
                            row = next_row
                            col = next_col
                            orbit_size += 1
                            break  # out of the for loop
                perm[row] = col
                logger.debug("set [" + str(row) + "] =" + str(col) + " (distance=" + str(record.distance) +")")
                used[col] = 1
                left -= 1
    return perm


def python_estimate_perm(op_type, op_order, molecule, dir):
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
            size = len(group)
            distances = np.zeros((size, size))
            for i in range(len(group)):
                for j in range(len(group)):
                    a = rotated[group[i]]
                    b = molecule.Q[group[j]]
                    distances[i][j] = np.linalg.norm(a - b)
        perm = perm_builder(op_type, op_order, group, distances, perm)
    return perm


def python_perm_builder(op_type, op_order, group, distances, perm):
    def recursive_perm_builder():

        # recursive:
        # find a pair of points (X,Y) minimally distant from each other (eg, minimum of distances_matrix
        # new_matrix = remove x row, y column from distances_matrix

        # INEFFICIENT METHOD:

        min_distance = MAXDOUBLE
        maxvals = np.ones(len(group)) * MAXDOUBLE
        left = len(group)
        while left > 0:
            if left < op_order:
                for index in left:
                    perm[index] = index
                    left -= 1
            else:
                flattened_index = np.argmin(
                    distances)  # numpy argmin returns the index of minimum value if matrix was flattened array
                i = int(flattened_index / len(group))  # this can be covnerted to an i,j index  by dividing and %ing
                j = int(flattened_index % len(group))
                perm[group[j]] = group[i]
                distances[:, j] = maxvals
                distances[i, :] = maxvals
                left -= 1

                # starting cycle
                orbitstart = i
                orbitdone = False
                orbitlength = 1
                while orbitlength < op_order:
                    (i, j) = (int(np.argmin(distances) / len(distances[0]))
                              , int(np.argmin(distances) % len(distances[0])))
                    perm[group[j]] = group[i]
                    distances[:, j] = maxvals
                    distances[i, :] = maxvals
                    left -= 1
                    orbitlength += 1
                perm[group[orbitstart]] = group[j]  # close loop
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


def test():
		for (j = 0; j < tableSize && left > 0; j++) {
			int enoughForFullOrbit = left >= opOrder;
			int row = distances[j].row;
			int col = distances[j].col;
			int dist;

			# If we have used this item already - skip it.
			if (perm[row] != -1)
				continue;

			# If we do not have enought to full groups, set all remaining items to themselves
			if (left == 1 || (type == CN && !enoughForFullOrbit)) {
				for (k = 0; k < groupSize; k++) {
					if (used[group[k]] == 0) {
						perm[group[k]] = group[k];
						LOG(debug) << "set [" << group[k] << "] =" << group[k] << " (distance=" << distances[k].distance << "???)";
					}
				}
				break;
			}
			if (opOrder == 2) {
				# Special treatment - only size 1 and two orbits are allowed
				# If both elements are not yet set, use the element.
				if (perm[row] == -1 && perm[col] == -1)  {
					perm[row] = col;
					perm[col] = row;
					LOG(debug) << "set [" << row << "] =" << col<< " (distance="<<distances[j].distance<<")";
					LOG(debug) << "set [" << col << "] =" << row << " (distance=" << distances[j].distance << "???)";
					left -= (row == col) ? 1 : 2;
				}
			}
			else {
				# we now want to complete an orbit.
				if (perm[row] == -1 && used[col] == 0) {
					perm[row] = col;
					used[col] = 1;
					LOG(debug) << "set [" << row << "] =" << col << " (distance=" << distances[j].distance << ")";
					left--;
				}
				else {
					continue;
				}

				# if this is an orbit of size one...
				if (row == col) continue;

				# If there is no more room for full orbit, must be SN and size two orbit
				if (type == SN && !enoughForFullOrbit) {
					perm[col] = row;
					used[row] = 1;
					LOG(debug) << "set [" << col << "] =" << row << " (distance=" << distances[j].distance << ")";
					left--;
					continue;
				}

				# Run until an orbit is complete
				orbitDone = false;
				orbitStart = row;
				orbitSize = 1;

				while (!orbitDone) {

					if (orbitSize == opOrder - 1) {
						row = col;
						col = orbitStart;
						orbitDone = true;
					}
					else {
						# Search for the next orbit element
						for (k = j + 1; k < tableSize; k++) {
							if (distances[k].row == col && used[distances[k].col] == 0 &&
								distances[k].col != distances[k].row) {
								if (orbitStart == distances[k].col) {
									if (type == SN && orbitSize == 1) {
										# we have now closed an orbit of size 2
										orbitDone = true;
									}
									else {
										continue;

									}
								}
								row = distances[k].row;
								col = distances[k].col;
								dist = distances[k].distance;

								orbitSize++;
								break;
							}
						}
					}
					perm[row] = col;
					used[col] = 1;
					LOG(debug) << "set [" << row << "] =" << col << " (distance=" << dist << ")";
					left--;
				}
			}
		}
		'''