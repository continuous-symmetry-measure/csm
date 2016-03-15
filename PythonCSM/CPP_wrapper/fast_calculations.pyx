import math

import numpy as np
cimport numpy as np
cimport csmlib
from calculations.constants import MAXDOUBLE, ZERO_IM_PART_MAX
from libcpp cimport bool

cdef class CalcState
cdef class Vector3D
cdef class Matrix3D

cdef build_polynomial(Vector3D lambdas, Vector3D m_t_B_2):  
    # The polynomial is described in equation 13.
    # The following code calculates the polynomial's coefficients quickly, and is taken
    # from the old C CSM code more or less untouched.
    cdef double coeffs[7]   # A polynomial of the 6th degree. coeffs[0] is for x^6, xoeefs[1] for x^5 , etc..
    coeffs[0]=1.0
    coeffs[1] = -2 * (lambdas.buf[0] + lambdas.buf[1] + lambdas.buf[2])
    coeffs[2] = lambdas.buf[0] * lambdas.buf[0] + lambdas.buf[1] * lambdas.buf[1] + lambdas.buf[2] * lambdas.buf[2] - \
                m_t_B_2.buf[0] - m_t_B_2.buf[1] - m_t_B_2.buf[2] + \
                4 * (lambdas.buf[0] * lambdas.buf[1] + lambdas.buf[0] * lambdas.buf[2] + lambdas.buf[1] * lambdas.buf[2])
    coeffs[3] = -8 * lambdas.buf[0] * lambdas.buf[1] * lambdas.buf[2] + \
                2 * (m_t_B_2.buf[0] * lambdas.buf[1] +
                     m_t_B_2.buf[0] * lambdas.buf[2] +
                     m_t_B_2.buf[1] * lambdas.buf[0] +
                     m_t_B_2.buf[1] * lambdas.buf[2] +
                     m_t_B_2.buf[2] * lambdas.buf[0] +
                     m_t_B_2.buf[2] * lambdas.buf[1] -
                     lambdas.buf[0] * lambdas.buf[2] * lambdas.buf[2] -
                     lambdas.buf[0] * lambdas.buf[0] * lambdas.buf[1] -
                     lambdas.buf[0] * lambdas.buf[0] * lambdas.buf[2] -
                     lambdas.buf[0] * lambdas.buf[1] * lambdas.buf[1] -
                     lambdas.buf[1] * lambdas.buf[1] * lambdas.buf[2] -
                     lambdas.buf[1] * lambdas.buf[2] * lambdas.buf[2])
    coeffs[4] = 4 * \
                ((lambdas.buf[0] * lambdas.buf[1] * lambdas.buf[2] * (lambdas.buf[0] + lambdas.buf[1] + lambdas.buf[2]) -
                  (m_t_B_2.buf[2] * lambdas.buf[0] * lambdas.buf[1] +
                   m_t_B_2.buf[1] * lambdas.buf[0] * lambdas.buf[2] +
                   m_t_B_2.buf[0] * lambdas.buf[2] * lambdas.buf[1]))) - \
                m_t_B_2.buf[0] * (lambdas.buf[1] * lambdas.buf[1] + lambdas.buf[2] * lambdas.buf[2]) - \
                m_t_B_2.buf[1] * (lambdas.buf[0] * lambdas.buf[0] + lambdas.buf[2] * lambdas.buf[2]) - \
                m_t_B_2.buf[2] * (lambdas.buf[0] * lambdas.buf[0] + lambdas.buf[1] * lambdas.buf[1]) + \
                lambdas.buf[0] * lambdas.buf[0] * lambdas.buf[1] * lambdas.buf[1] + \
                lambdas.buf[1] * lambdas.buf[1] * lambdas.buf[2] * lambdas.buf[2] + \
                lambdas.buf[0] * lambdas.buf[0] * lambdas.buf[2] * lambdas.buf[2]
    coeffs[5] = 2 * \
                (m_t_B_2.buf[0] * lambdas.buf[1] * lambdas.buf[2] * (lambdas.buf[1] + lambdas.buf[2]) +
                 m_t_B_2.buf[1] * lambdas.buf[0] * lambdas.buf[2] * (lambdas.buf[0] + lambdas.buf[2]) +
                 m_t_B_2.buf[2] * lambdas.buf[0] * lambdas.buf[1] * (lambdas.buf[0] + lambdas.buf[1])) \
                - 2 * \
                  (lambdas.buf[0] * lambdas.buf[1] * lambdas.buf[1] * lambdas.buf[2] * lambdas.buf[2] +
                   lambdas.buf[0] * lambdas.buf[0] * lambdas.buf[1] * lambdas.buf[2] * lambdas.buf[2] +
                   lambdas.buf[0] * lambdas.buf[0] * lambdas.buf[1] * lambdas.buf[1] * lambdas.buf[2])
    coeffs[6] = -m_t_B_2.buf[0] * lambdas.buf[1] * lambdas.buf[1] * lambdas.buf[2] * lambdas.buf[2] - \
                m_t_B_2.buf[1] * lambdas.buf[0] * lambdas.buf[0] * lambdas.buf[2] * lambdas.buf[2] - \
                m_t_B_2.buf[2] * lambdas.buf[0] * lambdas.buf[0] * lambdas.buf[1] * lambdas.buf[1] + \
                lambdas.buf[0] * lambdas.buf[0] * lambdas.buf[1] * lambdas.buf[1] * lambdas.buf[2] * lambdas.buf[2]

    return coeffs


def calculate_dir(bool is_zero_angle, int op_order, Vector3D lambdas, double lambda_max, Matrix3D m, Vector3D m_t_B, Vector3D B):
    cdef double m_max_B = 0.0
    cdef Vector3D dir = Vector3D.zero()
    cdef int i, j
    cdef double min_dist
    cdef int minarg

    # dir is calculated below according to formula (14) in the paper.
    # in the paper dir is called 'm_max'
    if is_zero_angle or op_order == 2:
        # If we are in zero teta case, we should pick the direction matching lambda_max
        min_dist = MAXDOUBLE
        minarg = 0

        for i in range(3):
            if math.fabs(lambdas.buf[i] - lambda_max) < min_dist:
                min_dist = math.fabs(lambdas.buf[i] - lambda_max)
                minarg = i
        for i in range(3):
            dir.buf[i] = m.buf[i][minarg]
    else:
        for i in range(3):
            for j in range(3):
                # error safety
                if math.fabs(lambdas.buf[j] - lambda_max) < 1e-6:
                    dir.buf[i] = m.buf[i][j]
                    break
                else:
                    dir.buf[i] += m_t_B.buf[j] / (lambdas.buf[j] - lambda_max) * m.buf[i][j]
            m_max_B = m_max_B + dir.buf[i] * B.buf[i]
    return dir, m_max_B


cdef PolynomialRoots(coeffs):
    cdef double coeffs_v[7]
    cdef double zeror[7]
    cdef double zeroi[7]

    cdef int i
    for i in range(7):
        coeffs_v[i] = coeffs[i]

    csmlib.rpoly(coeffs_v, 6, zeror, zeroi)
    cdef complex result[7]
    for i in range(7):
        result[i] = complex(zeror[i], zeroi[i])

    return result


cpdef get_lambda_max(Vector3D lambdas, Vector3D m_t_B_2):
    cdef double coeffs[7]
    coeffs = build_polynomial(lambdas, m_t_B_2)
    cdef complex roots[7]
    roots = PolynomialRoots(coeffs)
    # lambda_max is a real root of the polynomial equation
    # according to the description above the formula (13) in the paper
    cdef double lambda_max = -MAXDOUBLE
    cdef int i
    for i in range(len(roots)):
        if roots[i].real > lambda_max and math.fabs(roots[i].imag) < ZERO_IM_PART_MAX:
            lambda_max = roots[i].real
    return lambda_max

cpdef calc_ref_plane(int op_order, op_type, CalcState calc_state):
    cdef int i
    cdef Matrix3D m = Matrix3D()
    cdef Vector3D lambdas = Vector3D()
    csmlib.GetEigens(calc_state.A.buf, m.buf, lambdas.buf)
    cdef Vector3D m_t_B = m.T_mul_by_vec(calc_state.B)
    cdef Vector3D m_t_B_2 = Vector3D()
    for i in range(3):
        m_t_B_2[i] = m_t_B[i] * m_t_B[i]
    lambda_max=get_lambda_max(lambdas, m_t_B_2)
    dir, m_max_B = calculate_dir(op_type, op_order, lambdas, lambda_max, m, m_t_B, calc_state.B)
    csm = calc_state.CSM + (lambda_max - m_max_B) / 2
    csm = math.fabs(100 * (1.0 - csm / op_order))
    return csm, dir