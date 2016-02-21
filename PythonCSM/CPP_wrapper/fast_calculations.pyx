import math

import numpy as np
cimport numpy as np
cimport csmlib

from calculations.constants import MAXDOUBLE, ZERO_IM_PART_MAX

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

ITYPE = np.int
ctypedef np.int_t ITYPE_t

cimport cython
@cython.boundscheck(False)
def cross(np.ndarray[DTYPE_t, ndim=2] a, np.ndarray[DTYPE_t, ndim=2] b):
    cdef double *pa = <double *>a.data
    cdef double *pb = <double *>b.data
    return np.array([pa[1] * pb[2] - pa[2] * pb[1],
                     pa[2] * pb[0] - pa[0] * pb[2],
                     pa[0] * pb[1] - pa[1] * pb[0]]).T

def PolynomialRoots(coeffs):
    cdef double coeffs_v[7]
    cdef double zeror[7]
    cdef double zeroi[7]

    cdef int i
    for i in range(7):
        coeffs_v[i] = coeffs[i]

    csmlib.rpoly(coeffs_v, 6, zeror, zeroi)
    result = []
    for i in range(7):
        result.append(complex(zeror[i], zeroi[i]))

    return result


def calc_A_B(int op_order,
             int is_improper,
             np.ndarray[DTYPE_t, ndim=1] theta,
             np.ndarray[ITYPE_t, ndim=2] perms,
             int size,
             Q):
    # A is calculated according to formula (17) in the paper
    # B is calculated according to formula (12) in the paper

    A = np.zeros((3, 3,))
    B = np.zeros((1, 3,))  # Row vector for now

    # compute matrices according to current perm and its powers (the identity does not contribute anyway)
    cdef int i
    cdef int k

    for i in range(1, op_order):
        if is_improper and (i % 2):
            multiplier = -1 - math.cos(theta[i])
        else:
            multiplier = 1 - math.cos(theta[i])

        # The i'th power of the permutation
        cur_perm = perms[i]

        # Q_ is Q after applying the i'th permutation on atoms (Q' in the article)
        # Q_ = [Q[p] for p in cur_perm]  # Q'

        # A_intermediate is calculated according to the formula (5) in the paper
        for k in range(size):
            A = A + multiplier * ((Q[cur_perm[k]] @ Q[k].T) + (Q[k] @ Q[cur_perm[k]].T))
            B = B + math.sin(theta[i]) * cross(Q[k], Q[cur_perm[k]])

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


def calculate_csm(op_order, perms, size, Q, theta, lambda_max, m_max_B):
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
            dists += Q[k].T @ Q[cur_perm[k]]
        csm += math.cos(theta[i]) * dists

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
    theta = np.zeros(op_order)
    if not is_zero_angle:
        for i in range(1, op_order):
            theta[i] = 2 * math.pi * i / op_order

    perms = np.empty([op_order, size], dtype=np.int)
    perms[0] = [i for i in range(size)]
    for i in range(1, op_order):
        perms[i] = [perm[perms[i - 1][j]] for j in range(size)]

    # For all k, 0 <= k < size, Q[k] = column vector of x_k, y_k, z_k (position of the k'th atom)
    # - described on the first page of the paper
    Q = molecule.Q

    A, B = calc_A_B(op_order, is_improper, theta, perms, size, Q)

    # logger.debug("Computed matrix A is:")
    # logger.debug(A)
    # logger.debug("Computed vector B is: %s" % B)

    # lambdas - list of 3 eigenvalues of A
    # m - list of 3 eigenvectors of A
    lambdas, m = np.linalg.eig(A)
    # compute square of scalar multiplications of eigen vectors with B
    m_t_B = m.T @ B
    m_t_B_2 = np.power(m_t_B, 2)
    m_t_B_2 = m_t_B_2[:, 0]  # Convert from column vector to row vector

    # logger.debug("mTb: %s" % m_t_B)
    # logger.debug("mTb^2: %s" % m_t_B_2)



    coeffs = build_polynomial(lambdas, m_t_B_2)
    roots = np.roots(coeffs)
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

    csm = calculate_csm(op_order, perms, size, Q, theta, lambda_max, m_max_B)

    return csm, dir
