import math
import numpy

from permutations.lengths import len_molecule_permuter
from permutations.permuters import molecule_permuter
# from CPP_wrapper.fast_permutations import molecule_permuter
from calculations.molecule import ChainedPermutation
from CPP_wrapper import csm

__author__ = 'YAEL'

MAXDOUBLE = 100000000.0
ZERO_IM_PART_MAX = 1e-3


def perform_operation(csm_args, data):
    """
    Performs the csm calculation according to the options specified in csm_args
    :param csm_args: the csm arguments dictionary
    :param data: current calculations data object
    :return: the calculations data object with the permutation, csm, dMin and direction values updated
    """
    if 'dir' in csm_args:
        result = csm.FindBestPermUsingDir(data)
    else:
        if csm_args['findPerm']:
            result = csm.FindBestPerm(data)
        else:
            result = csm_operation(data, csm_args)  # csm_args['opName'], csm_args['molecule'].chains_perms)
    return result


def csm_operation(current_calc_data, csm_args):  # op_name, chains_perms):
    """
    Calculates minimal csm, dMin and directional cosines by applying permutations
    that keep the similar atoms within the group.
    Once it finds the optimal permutation , calls the CreateSymmetricStructure on the optimal permutation
    :param current_calc_data: current calculations data object
    :param args: The CSM arguments
    :return: the calculations data object with the permutation, csm, dMin and direction values updated
    """
    op_name = csm_args['opName']
    chained_perms = csm_args['molecule'].chained_perms

    result_csm = MAXDOUBLE
    dir = []
    optimal_perm = []
    current_calc_data.dir = [0, 0, 0]
    # calculate csm for each valid permutation & remember minimal

    if not chained_perms:
        # If no chained permutations specified - the regular permutations will be used
        chained_perms = [ChainedPermutation(1, list(range(len(current_calc_data.molecule.atoms))))]

    # Iterate through the permutations that swap chains
    for chained_perm in chained_perms:
        # and apply on each of them all the permutations on elements inside of each chain
        for perm in molecule_permuter(chained_perm.atom_perm,
                                      current_calc_data.molecule.equivalence_classes,
                                      current_calc_data.opOrder,
                                      current_calc_data.operationType == 'SN'):
            if csm_args['printPermutations']:
                print(perm)
            current_calc_data.perm = perm
            current_calc_data = calc_ref_plane(current_calc_data)
            # check, if it's a minimal csm, update dir and optimal perm
            if current_calc_data.csm < result_csm:
                (result_csm, dir, optimal_perm) = (current_calc_data.csm, current_calc_data.dir[:], perm[:])
            """print("CSM: %7lf\tPermutation: " % current_calc_data.csm + str(current_calc_data.perm)+"\tDirection: " +
                  str(current_calc_data.dir))"""

    if result_csm == MAXDOUBLE:
        # failed to find csm value for any permutation
        raise ValueError("Failed to calculate a csm value for %s" % op_name)

    current_calc_data.dMin = 1.0 - (result_csm / 100 * current_calc_data.opOrder / (current_calc_data.opOrder - 1))
    current_calc_data.perm = optimal_perm
    current_calc_data.dir = dir
    current_calc_data.csm = result_csm

    return csm.CreateSymmetricStructure(current_calc_data)


def total_number_of_permutations(csm_args):
    if csm_args['type'] != 'CH':
        return len_molecule_permuter(csm_args['molecule'], csm_args['opOrder'], csm_args['type'])
    else:
        num_perms = len_molecule_permuter(csm_args['molecule'], 2, 'SN')
        for i in range(2, csm_args['sn_max'] + 1, 2):
            num_perms += len_molecule_permuter(csm_args['molecule'], i, 'SN')
        return num_perms


def calc_ref_plane(current_calc_data):
    size = len(current_calc_data.molecule.atoms)
    is_improper = current_calc_data.operationType != 'CN'
    is_zero_angle = current_calc_data.operationType == 'CS'

    # TODO:	LOG(debug) << "calcRefPlane called";
    # TODO: LOG(debug) << "Permutation is " << permstrm.str();
    # TODO: LOG(debug) << "Direction is " << setprecision(2) << fixed << dir[0] << " " << dir[1] << " " << dir[2];

    print('calcRefPlane called')
    print('Permutation is ' + str(current_calc_data.perm))
    print("Direction is %lf %lf %lf" % (current_calc_data.dir[0], current_calc_data.dir[1], current_calc_data.dir[2]))

    # initialize identity permutation
    cur_perm = [i for i in range(size)]

    # For all k, 0 <= k < size, Q[k] = column vector of x_k, y_k, z_k (position of the k'th atom)
    # - described on the first page of the paper
    Q = [numpy.matrix([[atom.pos[0]], [atom.pos[1]], [atom.pos[2]]]) for atom in current_calc_data.molecule.atoms]
    Q_transpose = [Q[i].transpose() for i in range(size)]

    print("======================================================================\nQ:")
    for i in range(size):
        print("%d:" % i)
        print(Q[i])

    # A is calculated according to formula (17) in the paper
    A = numpy.matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])

    # B is calculated according to formula (12) in the paper
    B = numpy.matrix([0, 0, 0])

    numpy.set_printoptions(precision=4)


    # compute matrices according to current perm and its powers (the identity does not contribute anyway)
    for i in range(1, current_calc_data.opOrder):
        if is_zero_angle:
            theta = 0.0
        else:
            theta = 2 * math.pi * i / current_calc_data.opOrder

        # The i'th power of the permutation
        cur_perm = [current_calc_data.perm[cur_perm[j]] for j in range(size)]

        # Q_tag is Q after applying the i'th permutation on atoms
        Q_tag = [Q[cur_perm[i]] for i in range(size)]
        Q_transpose_tag = [Q_tag[i].transpose() for i in range(size)]

        print("Q_tag:")
        for j in range(size):
            print("%d:" % j)
            print(Q_tag[j])

        for j in range(size):
            print("%d:" % j)
            print(Q_transpose_tag[j])

        # A_intermediate is calculated according to the formula (5) in the paper
        A_intermediate = numpy.matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])

        for k in range(size):
            A_intermediate = A_intermediate + (Q_tag[k] * Q_transpose[k]) + (Q[k] * Q_transpose_tag[k])
            B = B + numpy.cross(Q[k].transpose(), Q_tag[k].transpose())

        B *= math.sin(theta)

        if is_improper and (i % 2) == 1:
            A = A + (A_intermediate * (-1 - math.cos(theta)))
        else:
            A = A + (A_intermediate * (1 - math.cos(theta)))

        print("Theta: %lf\tA - intermediate:" % theta)
        print(A_intermediate)

    # TODO: LOG(debug) << "Computed matrix is:" << setprecision(4);
    # TODO: LOG(debug) << A[0][0] << " " << A[0][1] << " " << A[0][2];
    # TODO: LOG(debug) << A[1][0] << " " << A[1][1] << " " << A[1][2];
    # TODO: LOG(debug) << A[2][0] << " " << A[2][1] << " " << A[2][2];
    print("Computed matrix is:")
    print(A)

    # lambdas - list of 3 eigenvalues of A
    # m - list of 3 eigenvectors of A
    lambdas, m = numpy.linalg.eig(A)
    # compute square of scalar multiplications of eigen vectors with B
    m_t_B = [(m[i].transpose() * B).tolist()[0][0] for i in range(3)]
    m_t_B_2 = [m_t_B[i] ** 2 for i in range(3)]

    # polynomial is calculated according to formula (13) in the paper

    # denominators[i] = (lambda_i - lambda_max)^2
    denominators = [numpy.polynomial.Polynomial([lambdas[i], -1]) ** 2 for i in range(3)]
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
    # TODO: LOG(debug) << "Coefficients: " << coeffs[0] << ", " << coeffs[1] << ", " << coeffs[2] << ", " << coeffs[3] << ", " << coeffs[4] << ", " << coeffs[5] << ", " << coeffs[6];
    print("Coefficients: ")
    print(polynomial)

    roots = polynomial.roots()

    # TODO: LOG(debug) << "rtr: " << roots[0] << " " << roots[1] << " " << roots[2] << " " << roots[3] << " " << roots[4]#  << " " << roots[5];
    print('roots: ')
    print(roots)

    # lambda_max is a real root of the polynomial equation
    # according to the description above the formula (13) in the paper
    lambda_max = -MAXDOUBLE
    for i in range(len(roots)):
        if roots[i].real > lambda_max and math.fabs(roots[i].imag) < ZERO_IM_PART_MAX:
            lambda_max = roots[i].real

    # TODO: LOG(debug) << setprecision(6) << fixed << "lambdas: " << lambdas[1] << " " << lambdas[2] << " " << lambdas[3];
    print("lambdas (eigenvalues): %lf %lf %lf" % (lambdas[0], lambdas[1], lambdas[2]))

    m_max_B = 0.0
    # dir is calculated below according to formula (14) in the paper.
    # in the paper dir is called 'm_max'
    if is_zero_angle or current_calc_data.opOrder == 2:
        # If we are in zero teta case, we should pick the direction matching lambda_max
        min_dist = MAXDOUBLE
        minarg = 0

        for i in range(3):
            if math.fabs(lambdas[i] - lambda_max) < min_dist:
                min_dist = math.fabs(lambdas[i] - lambda_max)
                minarg = i

        current_calc_data.dir = m[minarg].tolist()[0][:]

    else:
        for i in range(3):
            current_calc_data.dir[i] = 0.0
            for j in range(3):
                # error safety
                if math.fabs(lambdas[j] - lambda_max) < 1e-6:
                    current_calc_data.dir[i] = m[j][i]
                    break
                else:
                    current_calc_data.dir[i] += m_t_B[j] / (lambdas[j] - lambda_max) * m[j][i]
            m_max_B = m_max_B + current_calc_data.dir[i] * B[i]

    # initialize identity permutation
    cur_perm = [i for i in range(size)]

    # CSM is computed according to formula (15) in the paper with some changes according to CPP code:
    #  - lambda_max is divided by 2 (this appears only in CPP code)
    # - division by N (size) and by d^2 (mean square of Q - the original atom positions) is not
    #   made here but in the normalizations section
    current_calc_data.csm = 1.0
    for i in range(1, current_calc_data.opOrder):
        # This can be more efficient - save results of matrix computation
        if is_zero_angle:
            theta = 0.0
        else:
            theta = 2 * math.pi * i / current_calc_data.opOrder
        dists = 0.0

        # i'th power of permutation
        cur_perm = [current_calc_data.perm[cur_perm[j]] for j in range(size)]

        Q_tag = [Q[cur_perm[i]] for i in range(size)]

        for k in range(size):
            dists += (Q_transpose[k] * Q_tag[k]).tolist()[0][0]
        current_calc_data.csm += math.cos(theta) * dists

    # TODO: LOG(debug) << setprecision(6) << fixed << "csm=" << current_calc_data.csm << " lambda_max=" << lambda_max << " m_max_B=" << m_max_B;
    # TODO: LOG(debug) << setprecision(6) << fixed << "dir: " << current_calc_data.csm[0] << " " << current_calc_data.csm[1] << " " << current_calc_data.csm[2];
    print("csm=%lf lambda_max=%lf m_max_B=%lf" % (current_calc_data.csm, lambda_max, m_max_B))
    print("dir: %lf %lf %lf" % (current_calc_data.dir[0], current_calc_data.dir[1], current_calc_data.dir[2]))

    current_calc_data.csm += (lambda_max - m_max_B) / 2
    current_calc_data.csm = math.fabs(100 * (1.0 - current_calc_data.csm / current_calc_data.opOrder))

    # TODO: LOG(debug) << setprecision(6) << fixed << "dir - csm: " << current_calc_data.csm[0] << " " << current_calc_data.csm[1] << " " << current_calc_data.csm[2] << " - " << current_calc_data.csm;
    print("dir - csm: %lf %lf %lf - %lf" %
          (current_calc_data.dir[0], current_calc_data.dir[1], current_calc_data.dir[2], current_calc_data.csm))

    return current_calc_data
