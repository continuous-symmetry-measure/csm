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

class Operation:
    def __init__(self, op, sn_max=8):
        op=self._get_operation_data(op)
        self.type= op.type
        self.order = op.order
        if op.type=="CH":
            self.order=sn_max
        self.name = op.name

    def _get_operation_data(self, opcode):
        """
        Returns data about an operation based on the opcode
        Args:
            opcode: c2, s4, etc...

        Returns:
            And OperationCode object, with type, order and name
        """
        OperationCode = namedtuple('OperationCode', ('type', 'order', 'name'))
        _opcode_data = {
            "cs": ('CS', 2, "MIRROR SYMMETRY"),
            "ci": ('CI', 2, "INVERSION (S2)"),
            "ch": ('CH', 2, "CHIRALITY"),
            "c2": ('CN', 2, "C2 SYMMETRY"),
            'c3': ('CN', 3, "C3 SYMMETRY"),
            'c4': ('CN', 4, "C4 SYMMETRY"),
            'c5': ('CN', 5, "C5 SYMMETRY"),
            'c6': ('CN', 6, "C6 SYMMETRY"),
            'c7': ('CN', 7, "C7 SYMMETRY"),
            'c8': ('CN', 8, "C8 SYMMETRY"),
            'c10': ('CN', 10, "C10 SYMMETRY"),
            's1': ('CS', 2, "MIRROR SYMMETRY (S1)"),
            's2': ('SN', 2, "S2 SYMMETRY"),
            's4': ('SN', 4, "S4 SYMMETRY"),
            's6': ('SN', 6, "S6 SYMMETRY"),
            's8': ('SN', 8, "S8 SYMMETRY"),
            's10': ('SN', 8, "S10 SYMMETRY")
        }

        def isint(s):
            try:
                int(s)
                return True
            except ValueError:
                return False

        opcode = opcode.lower()
        if opcode[0] == 'c' and isint(opcode[1:]):
            return OperationCode(type='CN', order=int(opcode[1:]), name=opcode.upper() + ' SYMMETRY')
        if opcode[0] == 's' and isint(opcode[1:]):
            if opcode[1:] == '1':
                data = _opcode_data[opcode.lower()]
                return OperationCode(type=data[0], order=data[1], name=data[2])
            if int(opcode[1:]) % 2 != 0:
                raise ValueError("SN values must be even")
            return OperationCode(type='SN', order=int(opcode[1:]), name=opcode.upper() + ' SYMMETRY')
        try:
            data = _opcode_data[opcode.lower()]
        except KeyError:
            raise
        return OperationCode(type=data[0], order=data[1], name=data[2])

    @staticmethod
    def placeholder(op_type, op_order, sn_max=8):
        #make an arbitrary operation
        o=Operation("C2")
        #overwrite values to match input
        o.type=op_type
        o.order=op_order
        if op_type=="CH":
            o.order=sn_max
        return o

class Calculation:
    def __init__(self, operation, molecule, timeout=300, *args, **kwargs):
        self.operation=operation
        self.molecule=molecule
        self._timeout = timeout

    def calculate(self):
        pass

    @property
    def result(self):
        return self._csm_result

class CSMResult:
    """
    Holds the results of a CSM calculation. 
    When initialized, it handles some final processing of the result like denormalization and calculating local CSM.
    It is not intended to be initialized by an external user, but by the calculation.
    Properties:
    molecule: the original Molecule on which the calculation was called (denormalized)
    op_order: the order of the symmetry operation the calculation used
    op_type: the type of symmetry option (cs, cn, sn, ci, ch)
    csm: The final result CSM value
    perm: The final permutation that gave the result CSM value (a list of ints)
    dir: The final axis of symmetry
    perm_count: The number of permutations of the molecule
    local_csm: 
    d_min:
    symmetric_structure: the most symmetric structure when molecule permuted
    normalized_symmetric_structure: the normalized symmetric structure
    normalized_molecule_coords: the normalized molecule coordinates
    formula_csm: a recalculation of the symmetry value using a formula to measure symmetry. Should be close to identical
                    to algorithm result
    chain_perm: the chain permutation (by index-- names can be accessed with chain_perm_string)
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
        self._process_results()


    def _process_results(self):
        self.d_min = 1.0 - (self.csm / 100 * self.op_order / (self.op_order - 1))  # this is the scaling factor


        # save the normalized coords before we denormalize
        self.normalized_symmetric_structure = self.create_symmetric_structure(self.molecule, self.perm, self.dir, self.op_type,
                                                         self.op_order)
        self.normalized_molecule_coords=np.array(self.molecule.Q)

        #denormalize
        self.symmetric_structure = de_normalize_coords(list(self.normalized_symmetric_structure), self.molecule.norm_factor)
        self.molecule.de_normalize()

        # run and save the formula test
        self.formula_csm = self._formula_test()

        #get the chain perm
        self.get_chain_perm()

    def _formula_test(self):
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
        self.local_csm=local_csm
        return local_csm

    def get_chain_perm(self):
        '''
        finds the existing permutation between chains, in the result
        :return:
        '''
        molecule=self.molecule
        perm=self.perm
        chain_perm_dict = {}
        for chain in molecule.chains:
            index = molecule.chains[chain][0]
            permuted_index = perm[index]
            for chain2 in molecule.chains:
                if permuted_index in molecule.chains[chain2]:
                    chain_perm_dict[chain] = chain2
                    break
        chain_perm = []
        for chain in molecule.chains:
            permuted_index = chain_perm_dict[chain]
            chain_perm.append(permuted_index)
        self.chain_perm=chain_perm

    def chain_perm_string(self):
        chain_str=""
        for from_index, to_index in enumerate(self.chain_perm):
            from_chain=self.molecule.chains.index_to_string(from_index)
            to_chain=self.molecule.chains.index_to_string(to_index)
            chain_str+=from_chain + "->" + to_chain + ", "
        chain_str=chain_str[:-2] #remove final comma and space
        return chain_str



#TODO: replace all calls to this class with creation of Result class
def process_results(results, calc_local=False):
    """
    retained for legacy purposes
    """
    r=CSMResult(results, calc_local)
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