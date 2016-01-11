import csv

import math

import numpy as np
np.set_printoptions(precision=6)

# from permutations.lengths import len_molecule_permuter
from permutations.lengths import len_molecule_permuter
from permutations.permuters import MoleculePermuter, MoleculeLegalPermuter
from CPP_wrapper.fast_permutations import molecule_permuter
from calculations.molecule import ChainedPermutation
from CPP_wrapper import csm
import logging

logger = logging.getLogger("calculations")

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
    current_calc_data.dir = np.zeros(3)
    # calculate csm for each valid permutation & remember minimal

    csv_writer = None
    if csm_args['outPermFile']:
        csv_writer = csv.writer(csm_args['outPermFile'], lineterminator='\n')
        csv_writer.writerow(['Permutation', 'Direction', 'CSM'])

    if not chained_perms:
        # If no chained permutations specified - the regular permutations will be used
        chained_perms = [ChainedPermutation(1, list(range(len(current_calc_data.molecule.atoms))))]

    mp=MoleculeLegalPermuter(current_calc_data.molecule, current_calc_data.opOrder, current_calc_data.operationType == 'SN')
    # Iterate through the permutations that swap chains
    for chained_perm in chained_perms:
        # and apply on each of them all the permutations on elements inside of each chain
        for perm in mp.permute(): #chained_perm.atom_perm):
            current_calc_data.perm = perm
            current_calc_data = csm.CalcRefPlane(current_calc_data) # C++ version
            #current_calc_data = calc_ref_plane(current_calc_data) # Python version
            if current_calc_data.csm < result_csm:
                (result_csm, dir, optimal_perm) = (current_calc_data.csm, np.copy(current_calc_data.dir), perm[:])
            # check, if it's a minimal csm, update dir and optimal perm
            if csv_writer:
                csv_writer.writerow([perm, current_calc_data.dir, current_calc_data.csm])

    if result_csm == MAXDOUBLE:
        # failed to find csm value for any permutation
        raise ValueError("Failed to calculate a csm value for %s" % op_name)

    current_calc_data.dMin = 1.0 - (result_csm / 100 * current_calc_data.opOrder / (current_calc_data.opOrder - 1))
    current_calc_data.perm = optimal_perm
    current_calc_data.dir = dir
    current_calc_data.csm = result_csm

    #return csm.CreateSymmetricStructure(current_calc_data) #c++
    return create_symmetric_structure(current_calc_data) #python


def create_symmetric_structure(current_calc_data):
    def create_rotation_matrix(iOp):
        rot = np.zeros((3, 3))
        angle=0.0 if is_zero_angle else 2 * np.pi * iOp / current_calc_data.opOrder
        factor = -1 if is_improper and (iOp%2)==1 else 1

        # The rotation matrix is calculated similarly to the Rodrigues rotation matrix. The only
        # difference is that the matrix is also a reflection matrix when factor is -1.
        #
        # This is why we took the old C++ code instead of applying the Rodrigues formula directly.
        for s in range(3):
            for t in range(3):
                ang = np.cos(angle) if s==t else 0
                rot[s][t] = ang  + ((factor - np.cos(angle)) * m_dir[s] * m_dir[t] + np.sin(angle) * W[s][t])


        return rot

    logger.debug('***************************** Python ************************')
    logger.debug('createSymmetricStructure called')

    #set up some initial values:
    is_improper = current_calc_data.operationType != 'CN'
    is_zero_angle = current_calc_data.operationType == 'CS'

    cur_perm=np.arange(len(current_calc_data.perm)) #array of ints...

    m_pos=np.asarray([np.asarray(atom.pos) for atom in current_calc_data.molecule.atoms])
    current_calc_data.outAtoms=np.copy(m_pos)

    normalization=current_calc_data.dMin/current_calc_data.opOrder

    m_dir= current_calc_data.dir
    # The W matrix from the Rodrigues Rotation Formula (http://math.stackexchange.com/a/142831)
    W= np.array([[0.0, -m_dir[2], m_dir[1]], [m_dir[2], 0.0, -m_dir[0]], [-m_dir[1], m_dir[0], 0.0]])


    ########calculate and apply transform matrix#########
    ###for i<OpOrder
    for i in range(1,current_calc_data.opOrder):

        #get rotation
        rotation_matrix=create_rotation_matrix(i)
        rotated_positions = m_pos @ rotation_matrix

        #set permutation
        for w in range (len(cur_perm)):
            cur_perm[w]=current_calc_data.perm[cur_perm[w]]

        #add correct permuted rotation to atom in outAtoms
        for j in range (len(current_calc_data.outAtoms)):
            current_calc_data.outAtoms[j] += rotated_positions[cur_perm[j]]

    #apply normalization:
    current_calc_data.outAtoms *= normalization

    return current_calc_data






def total_number_of_permutations(csm_args):
    if csm_args['type'] != 'CH':
        return len_molecule_permuter(csm_args['molecule'], csm_args['opOrder'], csm_args['type'])
    else:
        num_perms = len_molecule_permuter(csm_args['molecule'], 2, 'SN')
        for i in range(2, csm_args['sn_max'] + 1, 2):
            num_perms += len_molecule_permuter(csm_args['molecule'], i, 'SN')
        return num_perms

def calc_refplane(current_calc_data):
    hi=0
    #arguments:
    #parray
    #perm
    #size
    #coef 3x3
    #multiplier

    #variables:
    #copyMat (3x3  mat)
    #"diag" (1X3 vec)
    #temp (1x3 vec)
    #"matrix" (3x3, all zeroes)
    #"vec" 1x3, all zeroes
    #cur perm- array of ints 12345..size
    #doubls csm, dists
    #bool isimporper
    #bool is zeroangle

    #step one: calc_A_B
    #permute permutation
    #call compute matrix
    #call compute vector

    #step two:
    #compute square of scalar multiplications of eigen vectors with b

    #step three:
    #build polynomial

    #step four:
    #solve polynomial

    #step five:
    #reset permutation to 123...size

    #step six:
    #compute CSM:


def calc_ref_plane(current_calc_data):
    size = len(current_calc_data.molecule.atoms)
    is_improper = current_calc_data.operationType != 'CN'
    is_zero_angle = current_calc_data.operationType == 'CS'

    logger.debug('***************************** Python ************************')
    logger.debug('calcRefPlane called')
    logger.debug('Permutation is ' + str(current_calc_data.perm))

    logger.debug("Direction is %lf %lf %lf" % (current_calc_data.dir[0], current_calc_data.dir[1], current_calc_data.dir[2]))

    # For all k, 0 <= k < size, Q[k] = column vector of x_k, y_k, z_k (position of the k'th atom)
    # - described on the first page of the paper
    def col_vec(list):
        a = np.array(list)
        a = a.reshape((3,1))
        return a

    Q = [col_vec(atom.pos) for atom in current_calc_data.molecule.atoms]

    def calc_A_B():
        # A is calculated according to formula (17) in the paper
        # B is calculated according to formula (12) in the paper

        cur_perm = [i for i in range(size)]
        A = np.zeros((3, 3,))
        B = np.zeros((1, 3)) # Row vector for now

        # compute matrices according to current perm and its powers (the identity does not contribute anyway)
        for i in range(1, current_calc_data.opOrder):
            if is_zero_angle:
                theta = 0.0
            else:
                theta = 2 * math.pi * i / current_calc_data.opOrder

            if is_improper and (i % 2):
                multiplier = -1 - math.cos(theta)
            else:
                multiplier = 1 - math.cos(theta)

            # The i'th power of the permutation
            cur_perm = [current_calc_data.perm[cur_perm[j]] for j in range(size)]


            # Q_ is Q after applying the i'th permutation on atoms (Q' in the article)
            Q_ = [Q[p] for p in cur_perm]  # Q'

            # A_intermediate is calculated according to the formula (5) in the paper
            for k in range(size):
                A = A + multiplier * ((Q_[k] @ Q[k].T) + (Q[k] @ Q_[k].T))
                B = B + math.sin(theta) * np.cross(Q[k].T, Q_[k].T)

        return A, B.T  # Return B as a column vector

    A,B = calc_A_B()

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
    if is_zero_angle or current_calc_data.opOrder == 2:
        # If we are in zero teta case, we should pick the direction matching lambda_max
        min_dist = MAXDOUBLE
        minarg = 0

        for i in range(3):
            if math.fabs(lambdas[i] - lambda_max) < min_dist:
                min_dist = math.fabs(lambdas[i] - lambda_max)
                minarg = i
        current_calc_data.dir = [m.tolist()[i][minarg] for i in range(3)]
    else:
        for i in range(3):
            current_calc_data.dir[i] = 0.0
            for j in range(3):
                # error safety
                if math.fabs(lambdas[j] - lambda_max) < 1e-6:
                    current_calc_data.dir[i] = m[i, j]
                    break
                else:
                    current_calc_data.dir[i] += m_t_B[j] / (lambdas[j] - lambda_max) * m[i, j]
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

        Q_ = [Q[cur_perm[i]] for i in range(size)]

        for k in range(size):
            dists += Q[k].T @ Q_[k]
        current_calc_data.csm += math.cos(theta) * dists

    logger.debug("csm=%lf lambda_max=%lf m_max_B=%lf" % (current_calc_data.csm, lambda_max, m_max_B))
    logger.debug("dir: %lf %lf %lf" % (current_calc_data.dir[0], current_calc_data.dir[1], current_calc_data.dir[2]))

    current_calc_data.csm += (lambda_max - m_max_B) / 2
    current_calc_data.csm = math.fabs(100 * (1.0 - current_calc_data.csm / current_calc_data.opOrder))

    logger.debug("dir - csm: %lf %lf %lf - %lf" %
          (current_calc_data.dir[0], current_calc_data.dir[1], current_calc_data.dir[2], current_calc_data.csm))

    return current_calc_data
