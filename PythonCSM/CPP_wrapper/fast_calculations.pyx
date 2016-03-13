import math

import numpy as np
cimport numpy as np
cimport csmlib
from calculations.constants import MAXDOUBLE, ZERO_IM_PART_MAX

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t


cpdef build_polynomial(lambdas, m_t_B_2):
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
