import csv
import math
import numpy as np
from calculations.constants import MINDOUBLE, MAXDOUBLE
from calculations.pair_cache import PairCache
from calculations.ref_plane import calc_ref_plane
from collections import namedtuple
from molecule.normalizations import de_normalize_coords, normalize_coords
from calculations.permuters import SinglePermPermuter, MoleculeLegalPermuter, CythonPermuter
import logging
from recordclass import recordclass

np.set_printoptions(precision=6)

logger = logging.getLogger("csm")

__author__ = 'YAEL'

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


def exact_calculation(op_type, op_order, molecule, perm=None, calc_local=False, permuter_class=MoleculeLegalPermuter, *args, **kwargs):
    if op_type == 'CH':  # Chirality
        sn_max = op_order
        # First CS
        best_result = csm_operation('CS', 2, molecule, perm, permuter_class)
        best_result.op_type='CS'
        if best_result.csm > MINDOUBLE:
            # Try the SN's
            for op_order in range(2, sn_max + 1, 2):
                result = csm_operation('SN', op_order, molecule, perm, permuter_class)
                if result.csm < best_result.csm:
                    best_result = result
                    best_result.op_type='SN'
                    best_result.op_order=op_order
                if best_result.csm > MINDOUBLE:
                    break

    else:
        best_result = csm_operation(op_type, op_order, molecule, perm, permuter_class)

    process_results(best_result, molecule)
    if calc_local:
        best_result.local_csm = compute_local_csm(molecule, best_result.perm, best_result.dir, best_result.op_type,
                                                  best_result.op_order)

    return best_result


def csm_operation(op_type, op_order, molecule, perm=None, permuter_class=MoleculeLegalPermuter):
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
        permuter = SinglePermPermuter(perm, molecule, op_order, op_type)
        logger.debug("SINGLE PERM")
    else:
        permuter = permuter_class(molecule, op_order, op_type)

    for pip in permuter.permute():
        csm, dir = calc_ref_plane(molecule, pip, op_order, op_type)
        if csm_state_tracer_func:
            traced_state.csm = csm
            traced_state.perm = pip.perm
            traced_state.dir = dir
            csm_state_tracer_func(traced_state)

        if csm < best_csm.csm:
            best_csm.csm = csm
            best_csm.dir = dir
            best_csm.perm = pip.perm[:]
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
