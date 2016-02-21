#
# Calculate the reference plane and CSM based on one permutation
#
# The calculations here are taken from the article Analytical Methods for Calculating Continuous Symmetry
# Measures and the Chirality Measure (Pinsky et al - 2008)
#
import logging
import numpy as np
from calculations.constants import ZERO_IM_PART_MAX, MAXDOUBLE
import math
# from CPP_wrapper.fast_calculations import cross
from CPP_wrapper import fast_calculations as cpp

from numpy.polynomial import Polynomial


# logger = logging.getLogger("csm")

def python_cross(a, b):
    return np.array([a[1] * b[2] - a[2] * b[1], a[2] * b[0]- a[0] * b[2],
                     a[0] * b[1] - a[1] * b[0]])
    return np.array([a[1][0] * b[2][0] - a[2][0] * b[1][0], a[2][0] * b[0][0] - a[0][0] * b[2][0],
                     a[0][0] * b[1][0] - a[1][0] * b[0][0]]).T

def python_outer_product(a,b):
    '''
    :param a: vector of length 3
    :param b: vector of length 3
    :return: 3x3 matrix of a@b's outer product
    '''
    return np.array([[a[0]*b[0], a[0]*b[1], a[0] * b[2]], [a[1]*b[0], a[1]*b[1], a[1] * b[2]], [a[2]*b[0], a[2]*b[1], a[2] * b[2]]])
    return np.array([[a[0][0]*b[0][0], a[0][0]*b[1][0], a[0][0] * b[2][0]], [a[1][0]*b[0][0], a[1][0]*b[1][0], a[1][0] * b[2][0]], [a[2][0]*b[0][0], a[2][0]*b[1][0], a[2][0] * b[2][0]]])

def python_inner_product(a,b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    return a[0][0]*b[0][0] + a[1][0]*b[1][0] + a[2][0]*b[2][0]

def python_matrix_by_vector(mat,v):
    x= mat[0][0]*v[0] + mat[0][1]*v[1] + mat[0][2]*v[2]
    y= mat[1][0]*v[0] + mat[1][1]*v[1] + mat[1][2]*v[2]
    z= mat[2][0]*v[0] + mat[2][1]*v[1] + mat[2][2]*v[2]
    #x= mat[0][0]*v[0][0] + mat[0][1]*v[1][0] + mat[0][2]*v[2][0]
    #y= mat[1][0]*v[0][0] + mat[1][1]*v[1][0] + mat[1][2]*v[2][0]
    #z= mat[2][0]*v[0][0] + mat[2][1]*v[1][0] + mat[2][2]*v[2][0]
    return np.array([x,y,z])

def calc_A_B(op_order, multiplier, sintheta, perms, size, Q):
    # A is calculated according to formula (17) in the paper
    # B is calculated according to formula (12) in the paper

    A = np.zeros((3, 3,))
    B = np.zeros((1, 3))  # Row vector for now

    # compute matrices according to current perm and its powers (the identity does not contribute anyway)
    for i in range(1, op_order):
        # The i'th power of the permutation
        cur_perm = perms[i]

        # Q_ is Q after applying the i'th permutation on atoms (Q' in the article)
        # Q_ = [Q[p] for p in cur_perm]  # Q'

        # A_intermediate is calculated according to the formula (5) in the paper
        for k in range(size):
            #A = A + multiplier[i] * ((Q[cur_perm[k]] @ Q[k].T) + (Q[k] @ Q[cur_perm[k]].T))
            A= A + multiplier[i] * (python_outer_product(Q[cur_perm[k]], Q[k]) + python_outer_product(Q[k], Q[cur_perm[k]]))
            B = B + sintheta[i] * np.cross(Q[k], Q[cur_perm[k]])

    return A, B.T  # Return B as a column vector


def build_polynomial(lambdas, m_t_B_2):
    # The polynomial is described in equation 13.
    # The following code calculates the polynomial's coefficients quickly, and is taken
    # from the old C CSM code more or less untouched.
    coeffs = [1.0, 0, 0, 0, 0, 0, 0]  # A polynomial of the 6th degree. coeffs[0] is for x^6, xoeefs[1] for x^5 , etc..
    coeffs[1] = -2 * (lambdas[0] + lambdas[1] + lambdas[2])
    coeffs[2] = lambdas[0] * lambdas[0] + lambdas[1] * lambdas[1] + lambdas[2] * lambdas[2] - \
                m_t_B_2[0] - m_t_B_2[1] - m_t_B_2[2] + \
                4 * (lambdas[0] * lambdas[1] + lambdas[0] * lambdas[2] + lambdas[1] * lambdas[2])
    coeffs[3] = -8 * lambdas[0] * lambdas[1] * lambdas[2] + \
                2 * (m_t_B_2[0] * lambdas[1] +
                     m_t_B_2[0] * lambdas[2] +
                     m_t_B_2[1] * lambdas[0] +
                     m_t_B_2[1] * lambdas[2] +
                     m_t_B_2[2] * lambdas[0] +
                     m_t_B_2[2] * lambdas[1] -
                     lambdas[0] * lambdas[2] * lambdas[2] -
                     lambdas[0] * lambdas[0] * lambdas[1] -
                     lambdas[0] * lambdas[0] * lambdas[2] -
                     lambdas[0] * lambdas[1] * lambdas[1] -
                     lambdas[1] * lambdas[1] * lambdas[2] -
                     lambdas[1] * lambdas[2] * lambdas[2])
    coeffs[4] = 4 * \
                ((lambdas[0] * lambdas[1] * lambdas[2] * (lambdas[0] + lambdas[1] + lambdas[2]) -
                  (m_t_B_2[2] * lambdas[0] * lambdas[1] +
                   m_t_B_2[1] * lambdas[0] * lambdas[2] +
                   m_t_B_2[0] * lambdas[2] * lambdas[1]))) - \
                m_t_B_2[0] * (lambdas[1] * lambdas[1] + lambdas[2] * lambdas[2]) - \
                m_t_B_2[1] * (lambdas[0] * lambdas[0] + lambdas[2] * lambdas[2]) - \
                m_t_B_2[2] * (lambdas[0] * lambdas[0] + lambdas[1] * lambdas[1]) + \
                lambdas[0] * lambdas[0] * lambdas[1] * lambdas[1] + \
                lambdas[1] * lambdas[1] * lambdas[2] * lambdas[2] + \
                lambdas[0] * lambdas[0] * lambdas[2] * lambdas[2]
    coeffs[5] = 2 * \
                (m_t_B_2[0] * lambdas[1] * lambdas[2] * (lambdas[1] + lambdas[2]) +
                 m_t_B_2[1] * lambdas[0] * lambdas[2] * (lambdas[0] + lambdas[2]) +
                 m_t_B_2[2] * lambdas[0] * lambdas[1] * (lambdas[0] + lambdas[1])) \
                - 2 * \
                  (lambdas[0] * lambdas[1] * lambdas[1] * lambdas[2] * lambdas[2] +
                   lambdas[0] * lambdas[0] * lambdas[1] * lambdas[2] * lambdas[2] +
                   lambdas[0] * lambdas[0] * lambdas[1] * lambdas[1] * lambdas[2])
    coeffs[6] = -m_t_B_2[0] * lambdas[1] * lambdas[1] * lambdas[2] * lambdas[2] - \
                m_t_B_2[1] * lambdas[0] * lambdas[0] * lambdas[2] * lambdas[2] - \
                m_t_B_2[2] * lambdas[0] * lambdas[0] * lambdas[1] * lambdas[1] + \
                lambdas[0] * lambdas[0] * lambdas[1] * lambdas[1] * lambdas[2] * lambdas[2]

    return coeffs


def calculate_dir(is_zero_angle, op_order, lambdas, lambda_max, m, m_t_B, B):
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
    return dir, m_max_B


def calculate_csm(op_order, perms, size, Q, costheta, lambda_max, m_max_B):
    # initialize identity permutation

    # CSM is computed according to formula (15) in the paper with some changes according to CPP code:
    #  - lambda_max is divided by 2 (this appears only in CPP code)
    # - division by N (size) and by d^2 (mean square of Q - the original atom positions) is not
    #   made here but in the normalizations section
    csm = 1.0
    for i in range(1, op_order):
        dists = 0.0

        # i'th power of permutation
        cur_perm = perms[i]

        # Q_ = [Q[cur_perm[i]] for i in range(size)]

        for k in range(size):
            dists += python_inner_product( Q[k], Q[cur_perm[k]])
        csm += costheta[i] * dists

    # logger.debug("csm=%lf lambda_max=%lf m_max_B=%lf" % (csm, lambda_max, m_max_B))
    # logger.debug("dir: %lf %lf %lf" % (dir[0], dir[1], dir[2]))

    csm += (lambda_max - m_max_B) / 2
    csm = math.fabs(100 * (1.0 - csm / op_order))

    # logger.debug("dir - csm: %lf %lf %lf - %lf" %
    #             (dir[0], dir[1], dir[2], csm))
    return csm


def calc_ref_plane(molecule, perm, op_order, op_type):
    size = len(molecule.atoms)
    is_improper = op_type != 'CN'
    is_zero_angle = op_type == 'CS'

    # pre-caching:
    multiplicand= 2 * math.pi /op_order
    costheta = np.zeros(op_order)
    sintheta = np.zeros(op_order)
    multiplier=np.zeros(op_order)
    if not is_zero_angle:
        for i in range(1, op_order):
            x= multiplicand * i
            costheta[i] = math.cos(x)
            sintheta[i] = math.sin(x)
            if is_improper and (i % 2):
                multiplier[i] = -1 - costheta[i]
            else:
                multiplier[i] = 1 - costheta[i]

    perms = np.empty([op_order, size], dtype=np.int)
    perms[0] = [i for i in range(size)]
    for i in range(1, op_order):
        perms[i] = [perm[perms[i - 1][j]] for j in range(size)]

    # For all k, 0 <= k < size, Q[k] = column vector of x_k, y_k, z_k (position of the k'th atom)
    # - described on the first page of the paper
    Q = molecule.Q

    A, B = calc_A_B(op_order, multiplier, sintheta, perms, size, Q)

    # logger.debug("Computed matrix A is:")
    # logger.debug(A)
    # logger.debug("Computed vector B is: %s" % B)

    # lambdas - list of 3 eigenvalues of A
    # m - list of 3 eigenvectors of A
    lambdas, m = np.linalg.eig(A)
    # compute square of scalar multiplications of eigen vectors with B
    #m_t_B = m.T @ B
    m_t_B=python_matrix_by_vector(m.T, B)
    m_t_B_2 = np.power(m_t_B, 2)
    #m_t_B_2 = m_t_B_2[:, 0]  # Convert from column vector to row vector

    # logger.debug("mTb: %s" % m_t_B)
    # logger.debug("mTb^2: %s" % m_t_B_2)



    coeffs = build_polynomial(lambdas, m_t_B_2)
    roots = cpp.PolynomialRoots(coeffs)
    # polynomial = build_polynomial()
    # roots = polynomial.roots()

    # logger.debug('roots: ')
    # logger.debug(roots)

    # lambda_max is a real root of the polynomial equation
    # according to the description above the formula (13) in the paper
    lambda_max = -MAXDOUBLE
    for i in range(len(roots)):
        if roots[i].real > lambda_max and math.fabs(roots[i].imag) < ZERO_IM_PART_MAX:
            lambda_max = roots[i].real

    # logger.debug("lambdas (eigenvalues): %lf %lf %lf" % (lambdas[0], lambdas[1], lambdas[2]))

    dir, m_max_B = calculate_dir(is_zero_angle, op_order, lambdas, lambda_max, m, m_t_B,B)

    csm = calculate_csm(op_order, perms, size, Q, costheta, lambda_max, m_max_B)

    return csm, dir
