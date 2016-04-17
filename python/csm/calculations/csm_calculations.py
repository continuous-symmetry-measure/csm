from datetime import datetime
import itertools

import math

import numpy as np
from csm.calculations.constants import MINDOUBLE, MAXDOUBLE
from csm.fast import calc_ref_plane

from collections import namedtuple
from csm.molecule.normalizations import de_normalize_coords, normalize_coords
from csm.fast import CythonPermuter, SinglePermPermuter, TruePermChecker, PQPermChecker, CythonPIP
import logging

np.set_printoptions(precision=6)

logger = logging.getLogger("csm")

__author__ = 'Itay, Devora, Yael'

CSMState = namedtuple('CSMState', ('molecule',
                                   'op_order',
                                   'op_type',
                                   'csm',
                                   'perm',
                                   'dir',
                                   'd_min',
                                   'symmetric_structure',
                                   'local_csm',
                                   'perm_count',))
CSMState.__new__.__defaults__ = (None,) * len(CSMState._fields)

# When this property is set by an outside caller, it is called every permutation iteration with the current CSMState
# This is useful for writing all permutations to file during the calculation
csm_state_tracer_func = None


def process_results(results):
    """
    Final normalizations and de-normalizations
    :param results: CSM old_calculations results
    """
    #    results.molecule.set_norm_factor(molecule.norm_factor)
    masses = [atom.mass for atom in results.molecule.atoms]
    normalize_coords(results.symmetric_structure, masses)

    results.molecule.de_normalize()
    symmetric_structure = de_normalize_coords(results.symmetric_structure, results.molecule.norm_factor)
    return results._replace(symmetric_structure=symmetric_structure)


def perm_count(op_type, op_order, molecule, keep_structure, print_perms=False, *args, **kwargs):
    def _len_group_permuter(group_size, cycle_size, add_cycles_of_two):
        """
        Returns the length of the group permuter.
        :param group_Size: Group size
        :param cycle_size: Cycle size
        :param add_cycles_of_two: When true, cycles of size of two are also legal
        :return: Number of permutations that will be returned by the group_permuter
        This is done by enumerating the different *kinds* of cycle_structs, and calculating the number of permutations
        for each kind.
        A *kind* of cycle struct is: (n1 cycles of size c1, n2 cycles of size c2, ...). Where c1..cn are either 1, cycle_size
        or 2 (if add_cycles_of_two is True), and the sum of ni*ci adds to group_size
        """

        def log_fact(n):
            """
            Returns the log(n!)
            """
            result = 0.0
            for i in range(1, n + 1):
                result += math.log(i)
            return result

        def log_nCr(n, r):
            """
            Calculates the log of n choose r
            """
            if r > n:
                raise ValueError("r must be <= n for nCr")
            if r < 0:
                raise ValueError("r must be non-negative")
            if n <= 0:
                raise ValueError("n must be positive")

            return log_fact(n) - log_fact(n - r) - log_fact(r)

        def generate_kinds():
            """
            Generates all cycle_struct kinds with the given group_size and cycle_lengths.
            A kind is a mapping from cycle_size to the number of cycles
            """
            add_two = add_cycles_of_two and cycle_size != 2  # Do we add cycles of two?
            for num_cycle_size in range(0, group_size // cycle_size + 1):
                if add_two:
                    left = group_size - cycle_size * num_cycle_size  # amount left
                    for num_two in range(0, left // 2 + 1):
                        yield {cycle_size: num_cycle_size, 2: num_two}
                else:
                    yield {cycle_size: num_cycle_size}

        def log_structs_of_kind(kind):
            """
            Returns the log of the number of cycle_structs of a given kind.
            :param kind: The kind of cycle_struct
            :return: The number of cycle_structs of kind.
            The calculation is done by choosing the cycles.

            For example, if there's only one cycle of size 2, there are group_size choose 2 different structs of this kind.
            If there are two cycles of size 2, there are group_size choose 2 * (group_size-2) choose 2 / 2
            The reasoning here is - first choose the first cycle, then choose the second cycle. Divide by two because the
            order of the two cycles doesn't matter.

            For 3 cycles of size 4 we get:
            group_size choose 4 * (group_size-4) choose 4 * (group_size-8) choose 4 / 3!
            We choose the first group of 4, then the second, then the third, but divide by 3! because the order of groups
            doesn't matter.

            If the group size is 12 and we allow 2 cycles of size 3 and 3 cycles of size 2, we get
            (12C3 * 9C3 / 2) * (6C2 * 4C2 * 2C2) / 2!
            We first choose the groups of size 3, then the groups of size 2.

            The above expressions can be further simplified, but we don't do it - this verison is more readable and performs
            very well.
            """
            log_result = 0
            size_left = group_size
            for cycle_size, cycle_count in kind.items():
                for i in range(cycle_count):
                    log_result += log_nCr(size_left, cycle_size)
                    size_left -= cycle_size
                log_result -= log_fact(cycle_count)
            return log_result

        def log_perms_in_kind(kind):
            """
            Returns the log of the number of permutations of a cycle_struct of a specific kind
            The number of permutations for a cycle of size n is (n-1)! so we just go and add that (add, because
            we're returning the log
            """
            log_result = 0
            for cycle_size, cycle_count in kind.items():
                log_result += cycle_count * log_fact(cycle_size - 1)
            return log_result

        count = 0
        for kind in generate_kinds():
            log_num_structs = log_structs_of_kind(kind)
            log_num_perms = log_perms_in_kind(kind)
            count += math.exp(log_num_structs + log_num_perms)
        return count

    perm_checker = TruePermChecker

    if not print_perms:
        if not keep_structure:
            count = 1
            groups = molecule.equivalence_classes
            for group in groups:
                count *= _len_group_permuter(len(group), op_order, op_type == 'SN')
            return int(count)

    if keep_structure:
        perm_checker = PQPermChecker

    traced_state = CSMState(molecule=molecule, op_type=op_type, op_order=op_order)
    permuter = CythonPermuter(molecule, op_order, op_type, perm_checker, perm_class=CythonPIP)

    for pip in permuter.permute():
        if csm_state_tracer_func:
            traced_state.csm = ''
            traced_state.perm = pip.perm
            traced_state.dir = ''
            csm_state_tracer_func(traced_state)
        if permuter.count % 1000000 == 0:
            print("counted", int(permuter.count / 1000000), "million permutations thus far...")
        count = permuter.count
    return count


def exact_calculation(op_type, op_order, molecule, keep_structure=False, perm=None, calc_local=False, *args, **kwargs):
    if keep_structure:
        permchecker = PQPermChecker
    else:
        permchecker = TruePermChecker

    permuter_class = CythonPermuter

    if op_type == 'CH':  # Chirality
        sn_max = op_order
        # First CS
        best_result = csm_operation('CS', 2, molecule, permuter_class, permchecker, perm)
        best_result.op_type = 'CS'
        if best_result.csm > MINDOUBLE:
            # Try the SN's
            for op_order in range(2, sn_max + 1, 2):
                result = csm_operation('SN', op_order, molecule, permuter_class, permchecker, perm)
                if result.csm < best_result.csm:
                    best_result = result
                    best_result.op_type = 'SN'
                    best_result.op_order = op_order
                if best_result.csm > MINDOUBLE:
                    break

    else:
        best_result = csm_operation(op_type, op_order, molecule, permuter_class, permchecker, perm)

    best_result = process_results(best_result)
    if calc_local:
        local_csm = compute_local_csm(molecule, best_result.perm, best_result.dir, best_result.op_type,
                                      best_result.op_order)
        best_result = best_result._replace(local_csm=local_csm)

    return best_result


def csm_operation(op_type, op_order, molecule, permuter_class=CythonPermuter, permchecker=TruePermChecker, perm=None):
    """
    Calculates minimal csm, dMin and directional cosines by applying permutations
    that keep the similar atoms within the group.
    Once it finds the optimal permutation , calls the CreateSymmetricStructure on the optimal permutation
    :param current_calc_data: current old_calculations data object
    :param args: The CSM arguments
    :return: A dictionary with all the results: csm, dMin, perm and direction
    """
    start_time = datetime.now()
    best_csm = CSMState(molecule=molecule, op_type=op_type, op_order=op_order, csm=MAXDOUBLE)
    traced_state = CSMState(molecule=molecule, op_type=op_type, op_order=op_order)

    if perm:
        permuter = SinglePermPermuter(np.array(perm), molecule, op_order, op_type)
        logger.debug("SINGLE PERM")
    else:
        permuter = permuter_class(molecule, op_order, op_type, permchecker)

    for calc_state in permuter.permute():
        if permuter.count % 1000000 == 0:
            print("calculated for", int(permuter.count / 1000000), "million permutations thus far...\t Time:",
                  datetime.now() - start_time)
        csm, dir = calc_ref_plane(op_order, op_type == 'CS', calc_state)
        if csm_state_tracer_func:
            traced_state.csm = csm
            traced_state.perm = calc_state.perm
            traced_state.dir = dir
            csm_state_tracer_func(traced_state)

        if csm < best_csm.csm:
            best_csm = best_csm._replace(csm=csm, dir=dir, perm=calc_state.perm)

    if best_csm.csm == MAXDOUBLE:
        # failed to find csm value for any permutation
        raise ValueError("Failed to calculate a csm value for %s" % op_type)

    best_csm = best_csm._replace(perm_count=permuter.count)
    d_min = 1.0 - (best_csm.csm / 100 * op_order / (op_order - 1))
    symmetric_structure = create_symmetric_structure(molecule, best_csm.perm, best_csm.dir, best_csm.op_type,
                                                     best_csm.op_order, d_min)
    best_csm = best_csm._replace(perm_count=permuter.count, d_min=d_min, symmetric_structure=symmetric_structure)
    return best_csm


def create_rotation_matrix(iOp, op_type, op_order, dir):
    is_improper = op_type != 'CN'
    is_zero_angle = op_type == 'CS'
    W = np.array([[0.0, -dir[2], dir[1]], [dir[2], 0.0, -dir[0]], [-dir[1], dir[0], 0.0]])
    rot = np.zeros((3, 3))
    angle = 0.0 if is_zero_angle else 2 * np.pi * iOp / op_order
    factor = -1 if is_improper and (iOp % 2) == 1 else 1

    # The rotation matrix is calculated similarly to the Rodrigues rotation matrix. The only
    # difference is that the matrix is also a reflection matrix when factor is -1.
    #
    # This is why we took the old C++ code instead of applying the Rodrigues formula directly.
    for s in range(3):
        for t in range(3):
            ang = np.cos(angle) if s == t else 0
            rot[s][t] = ang + ((factor - np.cos(angle)) * dir[s] * dir[t] + np.sin(angle) * W[s][t])

    return rot


def compute_local_csm(molecule, perm, dir, op_type, op_order):
    size = len(molecule.atoms)
    cur_perm = [i for i in range(size)]
    local_csm = np.zeros(size)
    m_pos = np.asarray([np.asarray(atom.pos) for atom in molecule.atoms])

    for i in range(op_order):
        rot = create_rotation_matrix(i, op_type, op_order, dir)

        # set permutation
        cur_perm = [perm[cur_perm[j]] for j in range(size)]

        # apply rotation to each atoms
        rotated = rot @ m_pos[cur_perm[i]]
        difference = rotated - m_pos[j]
        square = np.square(difference)
        sum = np.sum(square)
        local_csm[i] = sum * (100.0 / (2 * op_order))
    return local_csm


def create_symmetric_structure(molecule, perm, dir, op_type, op_order, d_min):
    logger.debug('create_symmetric_structure called')

    cur_perm = np.arange(len(perm))  # array of ints...
    size = len(perm)
    m_pos = np.asarray([np.asarray(atom.pos) for atom in molecule.atoms])
    symmetric = np.copy(m_pos)
    logger.debug("in atoms:")
    logger.debug(symmetric)

    normalization = d_min / op_order

    ########calculate and apply transform matrix#########
    ###for i<OpOrder
    for i in range(1, op_order):
        # get rotation
        rotation_matrix = create_rotation_matrix(i, op_type, op_order, dir)
        logger.debug("Rotation matrix:\n")
        logger.debug(rotation_matrix)
        # rotated_positions = m_pos @ rotation_matrix

        # set permutation
        cur_perm = [perm[cur_perm[j]] for j in range(size)]

        # add correct permuted rotation to atom in outAtoms
        for j in range(len(symmetric)):
            symmetric[j] += rotation_matrix @ m_pos[cur_perm[j]]
        logger.debug("Out atoms")
        logger.debug(symmetric)

    # apply normalization:
    symmetric *= normalization
    logger.debug("normalized out atoms:")
    logger.debug(symmetric)

    return symmetric


def approx_calculation(op_type, op_order, molecule, detect_outliers=True, *args, **kwargs):
    return find_best_perm(op_type, op_order, molecule, detect_outliers)


def find_best_perm(op_type, op_order, molecule, detect_outliers):
    logger.debug("findBestPerm called")

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
            curr_results = csm_operation(op_type, op_order, molecule, SinglePermPermuter, TruePermChecker, perm)
            best_for_this_dir = curr_results
            # (can add in conditions that curr_results.csm>1e-4 and curr_results.csm - old_resutls.csm>0.01 later)
            # (can also use max_iters instead of fixed 50, but 50 is whats given in the article and we don't care right now)
            for i in range(50):
                perm = estimate_perm(op_type, op_order, molecule, curr_results.dir)
                curr_results = csm_operation(op_type, op_order, molecule, SinglePermPermuter, TruePermChecker, perm)
                if curr_results.csm < best_for_this_dir.csm:
                    best_for_this_dir = curr_results

            if best_for_this_dir.csm < best.csm:
                best = best_for_this_dir

    return best


min_group_for_outliers = 10
def find_symmetry_directions(molecule, detect_outliers, op_type):
    # get average position of each equivalence group:
    group_averages = []
    for group in molecule.equivalence_classes:
        sum= [0,0,0]
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


def dir_fit(positions):
    def cross_product_minus_average(vector, average, matrix):
        '''
        :return: adds to matrix the cross_product of (vector-average) with itself.
        '''
        for i in range(3):
            for j in range(3):
                matrix[i][j] += vector[i] - average[i] * vector[j] - average[j]

    # get vector of the average position
    sum = np.einsum('ij->j', positions)
    average = sum / len(positions)

    mat = np.zeros((3, 3))
    for pos in positions:
        cross_product_minus_average(pos, average, mat)

        # computer eigenvalues and eigenvectors
        # lambdas - list of 3 eigenvalues of matrix
        # m - list of 3 eigenvectors of matrix
    lambdas, m = np.linalg.eig(mat)
    dirs = m

    # normalize result:
    for dir in dirs:
        normalize_dir(dir)

    return dirs


def normalize_dir(dir):
    dir=np.array(dir)
    norm = math.sqrt(math.fabs(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]))
    dir /= norm
    dir=list(dir)


def dirs_without_outliers(dirs, positions, op_type):
    additional_dirs=list(dirs)
    for dir in dirs:
        dists=[]
        # 1.Find the distance of each point from the line/plane
        for pos in positions:
            if op_type=='CS':
                dists.append(math.abs(dir[0]*pos[0]+dir[1]*pos[1]+dir[2]*pos[2]))
            else:
                dists.append(1)

        # 2. Find median of the distances m
        median=numpy.median(numpy.array(dists))
        # 3. for each distance di, if di > 2*m, remove it as an outlier
        with_outliers_removed=list()
        for i in range(len):
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
        elif math.fabs(dir[1] < MINDOUBLE):
            dir1 = [-dir[2], 0.0, dir[0]]
            dir2 = [0.0, 1.0, 0.0]
        else:
            dir1 = [-dir[1], dir[0], 0.0]
            dir2 = [0.0, -dir[2], dir[1]]
        # normalize dir1:
        normalize_dir(dir1)
        # remove projection of dir1 from dir2
        scal = dir1[0] * dir2[0] + dir1[1] * dir2[1] + dir1[2] * dir2[2]
        dir2[0] -= scal * dir1[0]
        dir2[1] -= scal * dir1[1]
        dir2[2] -= scal * dir1[2]
        # normalize dir2
        normalize_dir(dir2)
        added_dirs.append(dir1)
        added_dirs.append(dir2)

    return np.array(added_dirs)


def estimate_perm(op_type, op_order, molecule, dir):
    1



def blablablaperm():
    # step two: first estimation of the symmetry measure: first permutation
    # apply the symmetry operation once, using symm_element from step one

    # create rotation matrix
    rotation_mat = create_rotation_matrix(1, op_type, op_order, dir)

    # run rotation matrix on atoms
    rotated = (rotation_mat @ molecule.Q.T).T


    # create permutation:
    perm = [-1] * len(molecule)
    #permutation creation is done by group:
    for group in molecule.equivalence_classes:
        if len(group)==1:
            perm[group[0]]=group[0]

        else:
                # measure all the distances between all the points in X to all the points in Y
                distances = {}
                for i in group:
                    for j in group:
                        distances[(i,j)] = math.fabs(np.sum(rotated[i] - molecule.Q[j]))

                # recursive:
                # find a pair of points (X,Y) minimally distant from each other (eg, minimum of distances_matrix
                # new_matrix = remove x row, y column from distances_matrix

                #INEFFICIENT METHOD:

                while True:
                    min_distance = MAXDOUBLE
                    for (i,j) in distances:
                            if distances[(i,j)]<min_distance:
                                min_distance=distances[(i,j)]
                                pair=(i,j)
                    i=pair[0]
                    j=pair[1]
                    for k in group:
                        for t in group:
                            try:
                                if k==i:
                                    del distances[(k,t)]
                                elif t==j:
                                    del distances[(k, t)]
                            except:
                                pass

                    perm[i]=j
                    if min_distance==MAXDOUBLE:
                        break

    # return paired sets, as permutation
    return perm
