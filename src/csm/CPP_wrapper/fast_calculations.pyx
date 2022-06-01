#
# Calculate the reference plane and CSM based on one permutation
#
# The calculations here are taken from the article Analytical Methods for Calculating Continuous Symmetry
# Measures and the Chirality Measure (Pinsky et al - 2008)
#

from libc.math cimport fabs

import numpy as np
cimport numpy as np
cimport fastcpp
from csm.calculations.constants import MAXDOUBLE, ZERO_IM_PART_MAX
from libcpp cimport bool

cdef class CalcState
cdef class Vector3D
cdef class Matrix3D

np.set_printoptions(precision=20)

cdef build_polynomial(Vector3D lambdas, Vector3D m_t_B_2, double *coeffs):
    # The polynomial is described in equation 13.
    # The following code calculates the polynomial's coefficients quickly, and is taken
    # from the old C CSM code more or less untouched.
    # cdef double coeffs[7]   # A polynomial of the 6th degree. coeffs[0] is for x^6, xoeefs[1] for x^5 , etc..
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


def calculate_dir(bool is_zero_angle, int op_order, Vector3D lambdas, double lambda_max, Matrix3D m, Vector3D m_t_B, Vector3D B):
    cdef double m_max_B = 0.0
    cdef Vector3D dir = Vector3D.zero()
    cdef int i, j
    cdef double min_dist
    cdef int minarg

    # dir is calculated below according to formula (14) in the paper.
    # in the paper dir is called 'm_max'
    if is_zero_angle or op_order == 2:
        # If we are in zero theta case, we should pick the direction matching lambda_max
        min_dist = MAXDOUBLE
        minarg = 0

        for i in range(3):
            if fabs(lambdas.buf[i] - lambda_max) < min_dist:
                min_dist = fabs(lambdas.buf[i] - lambda_max)
                minarg = i
        for i in range(3):
            dir.buf[i] = m.buf[minarg][i]
    else:
        # print("Calculating direction")
        for i in range(3):
            dir.buf[i] = 0.0
            for j in range(3):
                # error safety
                if fabs(lambdas.buf[j] - lambda_max) < 1e-5:
                    dir.buf[i] = m.buf[j][i]
                    break
                else:
                    dir.buf[i] += m_t_B.buf[j] / (lambdas.buf[j] - lambda_max) * m.buf[j][i]
                #print("i=%d, j=%d" % (i, j))
                #print("dir[i] = %f" % dir.buf[i])

            #print("i=%d, dir[i] = %f" % (i, dir.buf[i]))
            m_max_B += dir.buf[i] * B.buf[i]

    #print("Returning direction ", dir[0], dir[1], dir[2])
    #print("Returning m_max_B ", m_max_B)
    return dir, m_max_B


cdef PolynomialRoots(double coeffs[7], complex *roots):
    cdef double zeror[6]
    cdef double zeroi[6]
    cdef int i

    fastcpp.rpoly(coeffs, 6, zeror, zeroi)
    for i in range(6):
        roots[i] = complex(zeror[i], zeroi[i])


cpdef get_lambda_max(Vector3D lambdas, Vector3D m_t_B_2, log=False):
    cdef double coeffs[7]
    #cdef double rounded_coeffs[7]
    cdef complex roots[6]
    cdef double lambda_max = -MAXDOUBLE
    cdef int i
    cdef int j
    log=False
    if log:
        print("get lambda max")
        print("lambdas", str(lambdas))
        print("m_t_B_2", str(m_t_B_2))

    if m_t_B_2[0] < ZERO_IM_PART_MAX and m_t_B_2[1] < ZERO_IM_PART_MAX and m_t_B_2[2] < ZERO_IM_PART_MAX:
        # In case m_t_B_2 is all zeros, we just get the maximum lambda
        lambda_max = lambdas[0]
        if lambdas[1] > lambda_max:
            lambda_max = lambdas[1]
        if lambdas[2] > lambda_max:
            lambda_max = lambdas[2]

        if log:
            print("m_t_B_2 is zero, returning the maximum lambda as is ", str(lambda_max))
        return lambda_max

    build_polynomial(lambdas, m_t_B_2, coeffs)

    roots=np.roots(coeffs)

    if log:
        print("coeffs")# and rounded")
        for i in range(7):
            print(coeffs[i]) #,"rounded:", rounded_coeffs[i])
        print("roots")
        for ro in roots:
            print(ro)

    # lambda_max is a real root of the polynomial equation
    # according to the description above the formula (13) in the paper
    for i in range(6):
        if roots[i].real > lambda_max and fabs(roots[i].imag) < ZERO_IM_PART_MAX:
            lambda_max = roots[i].real

    return lambda_max


def external_get_eigens(np.ndarray[DTYPE_t, ndim=2, mode="c"] A, np.ndarray[DTYPE_t, ndim=2, mode="c"] m, np.ndarray[DTYPE_t, ndim=1, mode="c"] lambdas):
   fastcpp.GetEigens( <double (*)[3]>A.data,  <double (*)[3]>m.data, <double *>lambdas.data)

cpdef calc_ref_plane(int op_order, bool is_op_cs, CalcState calc_state):
    global log
    cdef int i
    cdef int j


    #print("Perm:", str(calc_state.perms.get_perm(1)))
    log = False
    #if(list(calc_state.perm) ==[6,7,0,1,2,3,4,5]):
    #    log = True
    if log:
        print("Perm:")
        print(calc_state.perms.get_perm(1))
        print("A:")
        print(str(calc_state.A))

        print("B:")
        print(str(calc_state.B))

        print("preliminary CSM")
        print(str(calc_state.CSM))


    cdef Matrix3D m = Matrix3D.zero()
    cdef Vector3D lambdas = Vector3D.zero()
    fastcpp.GetEigens(calc_state.A.buf, m.buf, lambdas.buf)

    if log:
        print("m:")
        print(str(m))
        print("lambdas:")
        print(str(lambdas))

    cdef Vector3D m_t_B = m.T_mul_by_vec(calc_state.B)
    cdef Vector3D m_t_B_2 = Vector3D.zero()
    for i in range(3):
        m_t_B_2.buf[i] = m_t_B[i] * m_t_B[i]

    if log:
        print("m_t_B:")
        print(str(m_t_B))
        print("m_t_B_2:")
        print(str(m_t_B_2))

    lambda_max=get_lambda_max(lambdas, m_t_B_2, log)

    if log:
        print("lambda max")
        print(lambda_max)

    dir, m_max_B = calculate_dir(is_op_cs, op_order, lambdas, lambda_max, m, m_t_B, calc_state.B)

    if log:
        print("m_max_b:")
        print(m_max_B)
        print("dir:")
        print(str(dir))

    csm = calc_state.CSM + (lambda_max - m_max_B) / 2

    if log:
        print ("CSM step one (calc_state.CSM + (lambda_max - m_max_B) / 2)")
        print(str(csm))

    csm = fabs(100 * (1.0 - csm / op_order))

    if log:
        print ("CSM step two (fabs(100 * (1.0 - csm / op_order))")
        print("final-csm:")
        print(str(csm))

    return csm, dir.to_numpy()