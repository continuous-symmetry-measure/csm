import numpy as np
from csm.input_output.formatters import csm_log as print
from csm.calculations.basic_calculations import create_rotation_matrix, check_perm_cycles, \
    check_perm_structure_preservation
from csm.calculations.constants import MINDOUBLE, MAXDOUBLE
from csm.molecule.molecule import Molecule
from csm.molecule.normalizations import de_normalize_coords, normalize_coords
from collections import namedtuple
class CSMState(namedtuple('CSMState', ['molecule',
                                   'op_order',
                                   'op_type',
                                   'csm',
                                   'perm',
                                   'dir',
                                   'perm_count',
                                    'num_invalid'])):
    pass

CSMState.__new__.__defaults__ = (None,) * len(CSMState._fields)

class Operation:
    def __init__(self, op, sn_max=8, init=True):
        self.op_code=op
        if init:
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
        o=Operation("C2", init=False)
        #overwrite values to match input
        o.type=op_type
        o.order=op_order
        if op_type=="CH":
            o.order=sn_max
        return o

    def to_dict(self):
        return {
            "name":self.name,
            "order":self.order,
            "type":self.type
        }

    @staticmethod
    def from_dict(in_dict):
        #make an arbitrary operation
        o=Operation("C2", init=False)
        #overwrite values
        o.type=in_dict["type"]
        o.name=in_dict["name"]
        o.order=in_dict["order"]
        return o


class CSMResult:
    def __init__(self, state, operation, overall_stats={}, ongoing_stats={}):
        self.failed=False
        #input
        self.molecule=state.molecule.copy() #not yet denormalized
        self.normalized_molecule_coords = np.array(self.molecule.Q)
        self.molecule.de_normalize()
        self.operation=operation
        self.op_type=state.op_type
        self.op_order=state.op_order

        #result
        self.csm=state.csm
        self.dir=state.dir
        self.perm=state.perm
        self.normalized_symmetric_structure = self.create_symmetric_structure(self.normalized_molecule_coords, self.perm, self.dir, self.op_type,
                                                         self.op_order)
        self.symmetric_structure = de_normalize_coords(list(self.normalized_symmetric_structure), self.molecule.norm_factor)
        self.formula_csm = self.get_CSM_by_formula(self.molecule, self.symmetric_structure)

        #stats
        self.overall_statistics=overall_stats
        self.ongoing_statistics=ongoing_stats

        falsecount, num_invalid, cycle_counts, bad_indices = check_perm_cycles(self.perm, operation)
        self.overall_statistics["# bad cycles"]=falsecount
        self.overall_statistics["% bad cycles"]= num_invalid / len(self.molecule)
        try:
            self.overall_statistics["% structure"]=check_perm_structure_preservation(self.molecule, self.perm)
        except ValueError:
            self.overall_statistics["% structure"]= "n/a"

        if self.operation.name=="CHIRALITY":
            if self.op_type == 'CS':
                self.overall_statistics["best chirality"] = "CS"
            else:
                self.overall_statistics["best chirality"] =  "S%d" % (self.op_order)


        self.overall_statistics["formula CSM"]=self.formula_csm

    @property
    def d_min(self):
        return 1.0 - (self.csm / 100 * self.operation.order / (self.operation.order - 1))

    @property
    def local_csm(self):
        return self.compute_local_csm(self.molecule.Q, self.operation, self.dir)

    @property
    def chain_perm_string(self):
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
        chain_str=""
        for from_index, to_index in enumerate(chain_perm):
            from_chain=self.molecule.chains.index_to_string(from_index)
            to_chain=self.molecule.chains.index_to_string(to_index)
            chain_str+=from_chain + "->" + to_chain + ", "
        chain_str=chain_str[:-2] #remove final comma and space
        return chain_str

    def create_symmetric_structure(self, molecule_coords, perm, dir, op_type, op_order):
        # print('create_symmetric_structure called')

        cur_perm = np.arange(len(perm))  # array of ints...
        size = len(perm)
        m_pos = molecule_coords
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

    def get_CSM_by_formula(self, molecule, symmetric_structure):
        Q=molecule.Q
        # step one: get average of all atoms
        init_avg = np.mean(Q, axis=0)
        # step two: distance between intial and actual: initial - actual, squared
        # step three: normal: distance between initial and initial average, (x-x0)^2 + (y-y0)^2 + (z-z0)^2
        # step four: sum of distances between initial and actual, and then sum of x-y-z
        # step five: sum of normal
        distance = np.array([0.0, 0.0, 0.0])
        normal = 0.0
        for i in range(len(Q)):
            distance += (np.square(Q[i] - symmetric_structure[i]))  # square of difference
            normal += (np.sum(np.square(Q[i] - init_avg)))
        distance = np.sum(distance)
        # print("yaffa normal =", normal)
        # step six: 100 * step four / step five
        result = 100 * distance / normal
        return result

    def compute_local_csm(self, molecule_coords, operation, dir):
        size = len(molecule_coords)
        cur_perm = [i for i in range(size)]
        local_csm = np.zeros(size)
        m_pos = molecule_coords

        for i in range(operation.order):
            rot = create_rotation_matrix(i, operation.type, operation.order, dir)

            # set permutation
            cur_perm = [self.perm[cur_perm[j]] for j in range(size)]

            # apply rotation to each atoms
            rotated = rot @ m_pos[cur_perm[i]]
            difference = rotated - m_pos[i]
            square = np.square(difference)
            sum = np.sum(square)
            local_csm[i] = sum * (100.0 / (2 * operation.order))
        return local_csm

    def print_summary(self):
        try:
            percent_structure = check_perm_structure_preservation(self.molecule, self.perm)
            print("The permutation found maintains " +
                    str(round(percent_structure * 100, 2)) + "% of the original molecule's structure")

        except ValueError:
            print("The input molecule does not have bond information and therefore conservation of structure cannot be measured")

        falsecount, num_invalid, cycle_counts, bad_indices = check_perm_cycles(self.perm, self.operation)
        print(
            "The permutation found contains %d invalid %s. %.2lf%% of the molecule's atoms are in legal cycles" % (
                falsecount, "cycle" if falsecount == 1 else "cycles",
                    100 * (len(self.molecule) - num_invalid) / len(self.molecule)))

        for cycle_len in sorted(cycle_counts):
                valid = cycle_len == 1 or cycle_len == self.operation.order or (
                        cycle_len == 2 and self.operation.type == 'SN')
                count = cycle_counts[cycle_len]
                print("There %s %d %s %s of length %d" % (
                    "is" if count == 1 else "are", count, "invalid" if not valid else "",
                    "cycle" if count == 1 else "cycles",
                    cycle_len))

        if self.operation.name == "CHIRALITY":
            print("Minimum chirality was found in", self.overall_statistics["best chirality"])

        print("%s: %.4lf" % (self.operation.name, abs(self.csm)))

        print("Chain perm: " + self.chain_perm_string)


    def to_dict(self):
        return {"Result":
            {
                "molecule": self.molecule.to_dict(),
                "normalized_molecule_coords": [list(i) for i in self.normalized_molecule_coords],
                "operation": self.operation.to_dict(),

                "csm": self.csm,
                "perm": self.perm,
                "dir": list(self.dir),

                "normalized_symmetric_structure": [list(i) for i in self.normalized_symmetric_structure],
                "symmetric_structure": [list(i) for i in self.symmetric_structure],
                "formula_csm": self.formula_csm,

                "overall stats":self.overall_statistics,
                "ongoing stats":self.ongoing_statistics
            }
        }

    @staticmethod
    def from_dict():
        result_dict=input["Result"]
        molecule=Molecule.from_dict(result_dict["molecule"])
        molecule.normalize()
        operation=Operation.from_dict(result_dict["operation"])
        state=CSMState(molecule, operation.order, operation.type, result_dict["csm"], result_dict["perm"], result_dict["dir"])
        result=CSMResult(state, operation, result_dict["overall stats"], result_dict["ongoing stats"])
        return result

class FailedResult:
    def __init__(self, failed_reason, molecule, **kwargs):
        self.failed=True
        self.failed_reason=failed_reason

        self.molecule=molecule
        self.normalized_molecule_coords = []
        self.operation=kwargs["operation"]
        self.op_type=self.operation.type
        self.op_order=self.operation.order

        #result
        self.csm="n/a"
        self.dir=["n/a", "n/a", "n/a"]
        self.perm=["n/a"]
        self.normalized_symmetric_structure = []# [["n/a"] for i in range(len(molecule))]
        self.symmetric_structure = []#  [[0,0,0] for i in range(len(molecule))]
        self.formula_csm = "n/a"

        self.overall_statistics={
            "failed":"FAILED",
            "reason for failure":self.failed_reason
        }

        self.ongoing_statistics={}