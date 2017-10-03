import numpy as np
from csm.calculations.constants import MINDOUBLE
from csm.molecule.normalizations import de_normalize_coords, normalize_coords
from collections import namedtuple


class CSMState(namedtuple('CSMState', ['molecule',
                                   'op_order',
                                   'op_type',
                                   'csm',
                                   'perm',
                                   'dir',
                                   'perm_count'])):
    pass

CSMState.__new__.__defaults__ = (None,) * len(CSMState._fields)


class Result:
    """
    This class takes care of the final processing of the results, ie denormalization, calculating symmetric structure, etc
    """
    def __init__(self, state, calc_local=False):
        self.__CSMState=state
        self.molecule=state.molecule.copy()
        self.op_order=state.op_order
        self.op_type=state.op_type
        self.csm=state.csm
        self.perm=state.perm
        self.dir=state.dir
        self.perm_count=state.perm_count
        self.local_csm=""
        if calc_local:
            self.local_csm = self.compute_local_csm(self.molecule, self.perm, self.dir, self.op_type,
                                          self.op_order)
        #self.d_min=state.d_min
        #self.symmetric_structure=state.symmetric_structure
        #self.formula_csm=state.formula_csm
        #self.normalized_molecule_coords=state.normalized_molecule_coords
        #self.normalized_symmetric_structure=state.normalized_symmetric_structure
        self.process_results()


    def process_results(self):
        self.d_min = 1.0 - (self.csm / 100 * self.op_order / (self.op_order - 1))  # this is the scaling factor


        # save the normalized coords before we denormalize
        self.normalized_symmetric_structure = self.create_symmetric_structure(self.molecule, self.perm, self.dir, self.op_type,
                                                         self.op_order)
        self.normalized_molecule_coords=np.array(self.molecule.Q)

        #denormalize
        self.symmetric_structure = de_normalize_coords(self.normalized_symmetric_structure, self.molecule.norm_factor)
        self.molecule.de_normalize()

        # run and save the formula test
        self.formula_csm = self.formula_test()

    def formula_test(self):
        Q=self.molecule.Q
        # step one: get average of all atoms
        init_avg = np.mean(Q, axis=0)
        # step two: distance between intial and actual: initial - actual, squared
        # step three: normal: distance between initial and initial average, (x-x0)^2 + (y-y0)^2 + (z-z0)^2
        # step four: sum of distances between initial and actual, and then sum of x-y-z
        # step five: sum of normal
        distance = np.array([0.0, 0.0, 0.0])
        normal = 0.0
        for i in range(len(Q)):
            distance += (np.square(Q[i] - self.symmetric_structure[i]))  # square of difference
            normal += (np.sum(np.square(Q[i] - init_avg)))
        distance = np.sum(distance)
        # print("yaffa normal =", normal)
        # step six: 100 * step four / step five
        result = 100 * distance / normal
        return result

    def create_symmetric_structure(self, molecule, perm, dir, op_type, op_order):
        # print('create_symmetric_structure called')

        cur_perm = np.arange(len(perm))  # array of ints...
        size = len(perm)
        m_pos = np.asarray([np.asarray(atom.pos) for atom in molecule.atoms])
        symmetric = np.copy(m_pos)

        normalization = 1 / op_order

        ########calculate and apply transform matrix#########
        ###for i<OpOrder
        for i in range(1, op_order):
            # get rotation
            rotation_matrix = create_rotation_matrix(i, op_type, op_order, dir)
            # print("Rotation matrix:\n")
            # print(rotation_matrix)
            # rotated_positions = m_pos @ rotation_matrix

            # set permutation
            cur_perm = [perm[cur_perm[j]] for j in range(size)]

            # add correct permuted rotation to atom in outAtoms
            for j in range(len(symmetric)):
                symmetric[j] += rotation_matrix @ m_pos[cur_perm[j]]
                # print("Symmetric: ", symmetric)

        # apply normalization:
        symmetric *= normalization

        return symmetric

    def compute_local_csm(self, molecule, perm, dir, op_type, op_order):
        size = len(molecule)
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


#TODO: replace all calls to this class with creation of Result class
def process_results(results, calc_local=False):
    """
    retained for legacy purposes
    """
    r=Result(results, calc_local)
    return r


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






def check_perm_cycles(perm, op_order, op_type):
    '''
    This function checks the cycles in a given permutation according to the provided operation order and type.
    It counts the legal and illegal cycles in the permutation
    :param perm:
    :param op_order: 
    :param op_type: 
    :return: the number of illegal cycles, the number of molecules in illegal cycles, and a dictionary of cycle lengths, 
    with key =length cycle, val= mnumber of cycles of that length
    '''
    checked=[False]*len(perm)
    num_invalid=0
    truecount=0
    falsecount=0

    cycle_counts={}


    for i, index in enumerate(perm):
        if checked[i]:
            continue
        checked[i]=True
        cycle_len = 1
        while not checked[index]:
            checked[index] = True
            index = perm[index]
            cycle_len += 1

        if cycle_len in cycle_counts:
            cycle_counts[cycle_len]+=1
        else:
            cycle_counts[cycle_len]=1

        if cycle_len == 1 or cycle_len == op_order or (cycle_len == 2 and op_type == 'SN'):
            truecount += 1
        else:
            num_invalid += cycle_len
            falsecount+=1


    return falsecount, num_invalid, cycle_counts






def check_perm_equivalence(mol, perm):
    for origin, destination in enumerate(perm):
        if destination not in mol.atoms[origin].equivalency:
            return False
    return True

def check_perm_structure(mol, perm):
    if len(mol.bondset)==0:
        raise ValueError("Molecule does not have any bond information")

    broken=0
    for origin, destination in enumerate(perm):
        for adjacent in mol.atoms[origin].adjacent:
            if (destination, perm[adjacent]) not in mol.bondset:
                broken+=1

    percent_structure= (len(mol.bondset)- broken )/len(mol.bondset)

    return percent_structure


def array_distance(a, b):
    return np.sqrt(
        (a[0] - b[0]) * (a[0] - b[0])
        + (a[1] - b[1]) * (a[1] - b[1])
        + (a[2] - b[2]) * (a[2] - b[2]))