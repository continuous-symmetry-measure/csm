import numpy as np
from csm.calculations.constants import MINDOUBLE
from csm.molecule.normalizations import de_normalize_coords, normalize_coords
from collections import namedtuple

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


def process_results(results):
    """
    Final normalizations and de-normalizations
    :param results: CSM old_calculations results
    """
    #    results.molecule.set_norm_factor(molecule.norm_factor)
    d_min = 1.0 - (results.csm / 100 * results.op_order / (results.op_order - 1))

    symmetric_structure = create_symmetric_structure(results.molecule, results.perm, results.dir, results.op_type,
                                                     results.op_order, d_min)
    results = results._replace(d_min=d_min, symmetric_structure=symmetric_structure)
    masses = [atom.mass for atom in results.molecule.atoms]
    normalize_coords(results.symmetric_structure, masses)

    results.molecule.de_normalize()
    symmetric_structure = de_normalize_coords(results.symmetric_structure, results.molecule.norm_factor)
    results= results._replace(symmetric_structure=symmetric_structure)
    test2= yaffa_test(results)
    diff=results.csm - test2
    if diff>.0001 or diff<-.0001:
        print("LARGE DIFF:", diff, "yaffa", test2)
    else:
        print(diff)

    return results


def yaffa_test(result):
    #step one: get average of all atoms
    init_avg = np.mean(result.molecule.Q, axis=0)
    #step two: distance between intial and actual: initial - actual, squared
    #step three: normal: distance between initial and initial average, (x-x0)^2 + (y-y0)^2 + (z-z0)^2
    #step four: sum of distances between initial and actual, and then sum of x-y-z
    #step five: sum of normal
    distance=np.array([0.0, 0.0, 0.0])
    normal=0.0
    for i in range(len(result.molecule.Q)):
        initminusactual=result.molecule.Q[i]- result.symmetric_structure[i]
        square1=np.square(initminusactual)
        square2=initminusactual * initminusactual
        distance+=(np.square(result.molecule.Q[i]- result.symmetric_structure[i]))
        normal+=(np.sum(np.square(result.molecule.Q[i]-init_avg )))
    distance=np.sum(distance)
    #step six: 100 * step four / step five
    result= 100* distance /normal
    return result


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
        difference = rotated - m_pos[i]
        square = np.square(difference)
        sum = np.sum(square)
        local_csm[i] = sum * (100.0 / (2 * op_order))
    return local_csm


def create_symmetric_structure(molecule, perm, dir, op_type, op_order, d_min):
    #logger.debug('create_symmetric_structure called')

    cur_perm = np.arange(len(perm))  # array of ints...
    size = len(perm)
    m_pos = np.asarray([np.asarray(atom.pos) for atom in molecule.atoms])
    symmetric = np.copy(m_pos)

    normalization = d_min / op_order

    ########calculate and apply transform matrix#########
    ###for i<OpOrder
    for i in range(1, op_order):
        # get rotation
        rotation_matrix = create_rotation_matrix(i, op_type, op_order, dir)
     #   logger.debug("Rotation matrix:\n")
     #   logger.debug(rotation_matrix)
        # rotated_positions = m_pos @ rotation_matrix

        # set permutation
        cur_perm = [perm[cur_perm[j]] for j in range(size)]

        # add correct permuted rotation to atom in outAtoms
        for j in range(len(symmetric)):
            symmetric[j] += rotation_matrix @ m_pos[cur_perm[j]]

    # apply normalization:
    symmetric *= normalization

    return symmetric

