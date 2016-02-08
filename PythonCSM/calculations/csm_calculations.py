import csv
import math
import numpy as np
from collections import namedtuple
from molecule.normalizations import de_normalize_coords, normalize_coords
from calculations.permuters import MoleculePermuter, SinglePermPermuter
import logging
from recordclass import recordclass

np.set_printoptions(precision=6)

logger = logging.getLogger("csm")

__author__ = 'YAEL'

MINDOUBLE = 1e-8
MAXDOUBLE = 100000000.0
ZERO_IM_PART_MAX = 1e-3

CSMState = recordclass('CSMState', ('molecule',
                                    'op_order',
                                    'op_type',
                                    'csm',
                                    'perm',
                                    'dir',
                                    'd_min',
                                    'symmetric_structure',
                                    'local_csm'))
CSMState.__new__.__defaults__ = (None,) * len(CSMState._fields)

# When this property is set by an outside caller, it is called every permutation iteration with the current CSMState
# This is useful for writing all permutations to file during the calculation
csm_state_tracer_func = None

def process_results(results, keepCenter=False):
    """
    Final normalizations and de-normalizations
    :param results: CSM old_calculations results
    :param csm_args: CSM args
    """
    #    results.molecule.set_norm_factor(molecule.norm_factor)
    masses = [atom.mass for atom in results.molecule.atoms]
    normalize_coords(results.symmetric_structure, masses, keepCenter)

    results.molecule.de_normalize()
    results.symmetric_structure = de_normalize_coords(results.symmetric_structure, results.molecule.norm_factor)


def exact_calculation(op_type, op_order, molecule, perm=None, calc_local=False, permuter_class=MoleculePermuter, *args, **kwargs):
    if op_type == 'CH':  # Chirality
        sn_max = op_order
        # First CS
        best_result = csm_operation('CS', 2, molecule, perm, permuter_class)
        if best_result.csm > MINDOUBLE:
            # Try the SN's
            for op_order in range(2, sn_max + 1, 2):
                result = csm_operation('SN', op_order, molecule, perm, permuter_class)
                if result.csm < best_result.csm:
                    best_result = result
            
    else:
        best_result = csm_operation(op_type, op_order, molecule, perm, permuter_class)

    process_results(best_result, molecule)
    if calc_local:
        best_result.local_csm = compute_local_csm(molecule, best_result.perm, best_result.dir, best_result.op_type,
                                                  best_result.op_order)

    return best_result


def csm_operation(op_type, op_order, molecule, perm=None, permuter_class=MoleculePermuter):
    """
    Calculates minimal csm, dMin and directional cosines by applying permutations
    that keep the similar atoms within the group.
    Once it finds the optimal permutation , calls the CreateSymmetricStructure on the optimal permutation
    :param current_calc_data: current old_calculations data object
    :param args: The CSM arguments
    :return: A dictionary with all the results: csm, dMin, perm and direction
    """
    logger.debug("csm_op atoms:")
    logger.debug([atom.pos for atom in molecule.atoms])
    best_csm = CSMState(molecule=molecule, op_type=op_type, op_order=op_order, csm=MAXDOUBLE)
    traced_state = CSMState(molecule=molecule, op_type=op_type, op_order=op_order)

    if perm:
        permuter = SinglePermPermuter(perm)
        logger.debug("SINGLE PERM")
    else:
        permuter = permuter_class(molecule, op_order, op_type == 'SN')

    for perm in permuter.permute():
        csm, dir = calc_ref_plane(molecule, perm, op_order, op_type)
        if csm_state_tracer_func:
            traced_state.csm = csm
            traced_state.perm = perm
            traced_state.dir = dir
            csm_state_tracer_func(traced_state)

        if csm < best_csm.csm:
            best_csm.csm = csm
            best_csm.dir = dir
            best_csm.perm = perm[:]
            # TODO: Write permutations while looping

    if best_csm.csm == MAXDOUBLE:
        # failed to find csm value for any permutation
        raise ValueError("Failed to calculate a csm value for %s" % op_type)

    best_csm.d_min = 1.0 - (best_csm.csm / 100 * op_order / (op_order - 1))

    best_csm.symmetric_structure = create_symmetric_structure(molecule, best_csm.perm, best_csm.dir, best_csm.op_type,
                                                              best_csm.op_order, best_csm.d_min)
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


def calc_ref_plane(molecule, perm, op_order, op_type):
    size = len(molecule.atoms)
    is_improper = op_type != 'CN'
    is_zero_angle = op_type == 'CS'

    logger.debug('************* calc_ref_plane Python ************************')
    logger.debug('Permutation is ' + str(perm))

#    logger.debug("calcrefplane atoms:")
#    logger.debug([atom.pos for atom in molecule.atoms])

    # For all k, 0 <= k < size, Q[k] = column vector of x_k, y_k, z_k (position of the k'th atom)
    # - described on the first page of the paper
    def col_vec(list):
        a = np.array(list)
        a = a.reshape((3, 1))
        return a

    Q = [col_vec(atom.pos) for atom in molecule.atoms]

    def calc_A_B():
        # A is calculated according to formula (17) in the paper
        # B is calculated according to formula (12) in the paper

        cur_perm = [i for i in range(size)]
        A = np.zeros((3, 3,))
        B = np.zeros((1, 3))  # Row vector for now

        # compute matrices according to current perm and its powers (the identity does not contribute anyway)
        for i in range(1, op_order):
            if is_zero_angle:
                theta = 0.0
            else:
                theta = 2 * math.pi * i / op_order

            if is_improper and (i % 2):
                multiplier = -1 - math.cos(theta)
            else:
                multiplier = 1 - math.cos(theta)

            # The i'th power of the permutation
            cur_perm = [perm[cur_perm[j]] for j in range(size)]

            # Q_ is Q after applying the i'th permutation on atoms (Q' in the article)
            Q_ = [Q[p] for p in cur_perm]  # Q'

            # A_intermediate is calculated according to the formula (5) in the paper
            for k in range(size):
                A = A + multiplier * ((Q_[k] @ Q[k].T) + (Q[k] @ Q_[k].T))
                B = B + math.sin(theta) * np.cross(Q[k].T, Q_[k].T)

        return A, B.T  # Return B as a column vector

    A, B = calc_A_B()

    logger.debug("Computed matrix A is:")
    logger.debug(A)
    logger.debug("Computed vector B is: %s" % B)

    # lambdas - list of 3 eigenvalues of A
    # m - list of 3 eigenvectors of A
    lambdas, m = np.linalg.eig(A)
    # compute square of scalar multiplications of eigen vectors with B
    m_t_B = m.T @ B
    m_t_B_2 = np.power(m_t_B, 2)

    logger.debug("mTb: %s" % m_t_B)
    logger.debug("mTb^2: %s" % m_t_B_2)

    # polynomial is calculated according to formula (13) in the paper

    # denominators[i] = (lambda_i - lambda_max)^2
    denominators = [np.polynomial.Polynomial([lambdas[i], -1]) ** 2 for i in range(3)]
    numerator = m_t_B_2[0] * denominators[1] * denominators[2] + \
                m_t_B_2[1] * denominators[0] * denominators[2] + \
                m_t_B_2[2] * denominators[0] * denominators[1]
    denominator = denominators[0] * denominators[1] * denominators[2]
    if denominator != 0:
        # the polynomial is "numerator / denominator = 1"
        # multiply both sides by denominator and get: "numerator = denominator"
        # iff "denominator - numerator = 0"
        polynomial = denominator - numerator
    else:
        raise ValueError("Can't compute the polynomial - division by zero")

    # solve polynomial and find maximum eigenvalue and eigen vector
    logger.debug("Coefficients: ")
    logger.debug(polynomial)

    roots = polynomial.roots()

    logger.debug('roots: ')
    logger.debug(roots)

    # lambda_max is a real root of the polynomial equation
    # according to the description above the formula (13) in the paper
    lambda_max = -MAXDOUBLE
    for i in range(len(roots)):
        if roots[i].real > lambda_max and math.fabs(roots[i].imag) < ZERO_IM_PART_MAX:
            lambda_max = roots[i].real

    logger.debug("lambdas (eigenvalues): %lf %lf %lf" % (lambdas[0], lambdas[1], lambdas[2]))

    m_max_B = 0.0
    # dir is calculated below according to formula (14) in the paper.
    # in the paper dir is called 'm_max'
    if is_zero_angle or op_order == 2:
        # If we are in zero teta case, we should pick the direction matching lambda_max
        min_dist = MAXDOUBLE
        minarg = 0

        for i in range(3):
            if math.fabs(lambdas[i] - lambda_max) < min_dist:
                min_dist = math.fabs(lambdas[i] - lambda_max)
                minarg = i
        dir = [m.tolist()[i][minarg] for i in range(3)]
    else:
        dir = np.zeros(3)
        for i in range(3):
            for j in range(3):
                # error safety
                if math.fabs(lambdas[j] - lambda_max) < 1e-6:
                    dir[i] = m[i, j]
                    break
                else:
                    dir[i] += m_t_B[j] / (lambdas[j] - lambda_max) * m[i, j]
            m_max_B = m_max_B + dir[i] * B[i]

    # initialize identity permutation
    cur_perm = [i for i in range(size)]

    # CSM is computed according to formula (15) in the paper with some changes according to CPP code:
    #  - lambda_max is divided by 2 (this appears only in CPP code)
    # - division by N (size) and by d^2 (mean square of Q - the original atom positions) is not
    #   made here but in the normalizations section
    csm = 1.0
    for i in range(1, op_order):
        # This can be more efficient - save results of matrix computation
        if is_zero_angle:
            theta = 0.0
        else:
            theta = 2 * math.pi * i / op_order
        dists = 0.0

        # i'th power of permutation
        cur_perm = [perm[cur_perm[j]] for j in range(size)]

        Q_ = [Q[cur_perm[i]] for i in range(size)]

        for k in range(size):
            dists += Q[k].T @ Q_[k]
        csm += math.cos(theta) * dists

    logger.debug("csm=%lf lambda_max=%lf m_max_B=%lf" % (csm, lambda_max, m_max_B))
    logger.debug("dir: %lf %lf %lf" % (dir[0], dir[1], dir[2]))

    csm += (lambda_max - m_max_B) / 2
    csm = math.fabs(100 * (1.0 - csm / op_order))

    logger.debug("dir - csm: %lf %lf %lf - %lf" %
                 (dir[0], dir[1], dir[2], csm))

    return csm, dir
