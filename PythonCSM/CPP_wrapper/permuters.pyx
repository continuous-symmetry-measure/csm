import itertools

import math

import numpy as np
cimport numpy as np
from calculations.pair_cache import PairCache

__author__ = 'Devora'


class PermChecker:
    def __init__(self, mol):
        self.mol = mol

    def is_legal(self, pip, origin, destination):
        for adjacent in self.mol.atoms[destination].adjacent:
            if pip.p[adjacent] != -1 and (origin, pip.p[adjacent]) not in self.mol.bondset:
                return False
        return True


class PQPermChecker:
    def __init__(self, mol):
        self.mol = mol

    def is_legal(self, pip, origin, destination):
        for adjacent in self.mol.atoms[destination].adjacent:
            if pip.p[adjacent] != -1 and (origin, pip.p[adjacent]) not in self.mol.bondset:
                return False
            for adjacent in self.mol.atoms[origin].adjacent:
                if pip.q[adjacent] != -1 and (destination, pip.q[adjacent]) not in self.mol.bondset:
                    return False
        return True


class TruePermChecker:
    def __init__(self, mol):
        pass

    def is_legal(self, pip, origin, destination):
        return True

cdef _partial_calculate(self, group, cache):
        '''
        :param group: the group that was just permuted. represents the indexes in self.perm that need to have A,B calculated
        '''
        # permuted_array=array[perm[i]]
        # perms[i] = [perm[perms[i - 1][j]] for j in range(size)]
        cdef int iop
        cdef int j
        cdef int index
        cdef int permuted_index
        for iop in range(1, self.op_order):
            dists=0.0
            for j in range(len(group)):
                index = group[j]
                permuted_index=self._calc.perms[iop - 1][self.p[index]]
                self._calc.perms[iop][index] = permuted_index
                self._calc.A+=self.multiplier[iop] * cache.outer_product_sum(index, permuted_index)
                self._calc.B+=self.sintheta[iop]*cache.cross(index, permuted_index)
                dists += cache.inner_product(index, permuted_index)
            self._calc.CSM += self.costheta[iop] * dists

def partial_calculate(self, group, cache):
        _partial_calculate(self, group, cache)

class PartialCalculation:
        def __init__(self,A,B,CSM,perms):
            self.A=np.copy(A)
            self.B=np.copy(B)
            self.CSM=CSM
            self.perms=np.copy(perms)

        @classmethod
        def copyconstruct(cls,partcalc):
            return cls(partcalc.A,partcalc.B,partcalc.CSM,partcalc.perms)

        @classmethod
        def initialConstructor(cls, op_order, size):
            perms = -1 * np.ones([op_order, size], dtype=np.int)
            perms[0] = [i for i in range(size)]  # identity perm
            A = np.zeros((3, 3,))
            B = np.zeros((3,), dtype=np.float64, order="c")
            CSM=1.0
            return cls(A,B,CSM,perms)

class CABPermInProgress:
    class PartialCalculation:
        def __init__(self,A,B,CSM,perms):
            self.A=np.copy(A)
            self.B=np.copy(B)
            self.CSM=CSM
            self.perms=np.copy(perms)

        @classmethod
        def copyconstruct(cls,partcalc):
            return cls(partcalc.A,partcalc.B,partcalc.CSM,partcalc.perms)

        @classmethod
        def initialConstructor(cls, op_order, size):
            perms = -1 * np.ones([op_order, size], dtype=np.int)
            perms[0] = [i for i in range(size)]  # identity perm
            A = np.zeros((3, 3,))
            B = np.zeros((3,), dtype=np.float64, order="c")
            CSM=1.0
            return cls(A,B,CSM,perms)

    def __init__(self, mol, op_order, op_type, permchecker=PQPermChecker):
        size = len(mol.atoms)
        self.type="AB_cython"
        self.size=size
        self.p = [-1] * size
        self.q = [-1] * size
        self.permchecker = permchecker(mol)
        self.mol = mol
        self.op_order = op_order
        self.sintheta, self.costheta, self.multiplier, self.is_zero_angle = self.precalculate(op_type, op_order)
        self._calc=self.PartialCalculation.initialConstructor(op_order,size)

    @property
    def perm(self):
        return self.p

    @property
    def B(self):
        return self._calc.B

    @property
    def A(self):
        return self._calc.A

    @property
    def perms(self):
        return self._calc.perms

    @property
    def CSM(self):
        return self._calc.CSM

    def switch(self, origin, destination):
        if self.permchecker.is_legal(self, origin, destination):
            assert self.p[origin] == -1 and self.q[destination] == -1
            self.p[origin] = destination
            self.q[destination] = origin
            return True
        return False

    def unswitch(self, origin, destination):
        assert self.p[origin] == destination and self.q[destination] == origin
        self.p[origin] = -1
        self.q[destination] = -1

    def precalculate(self, op_type, op_order):
        is_improper = op_type != 'CN'
        is_zero_angle = op_type == 'CS'
        # pre-caching:
        sintheta = np.zeros(op_order)
        costheta = np.zeros(op_order)
        multiplier = np.zeros(op_order)

        for i in range(1, op_order):
            if not is_zero_angle:
                theta = 2 * math.pi * i / op_order
            cos=math.cos(theta)
            costheta[i]=cos
            sintheta[i]=math.sin(theta)
            if is_improper and (i % 2):
                multiplier[i] = -1 - cos
            else:
                multiplier[i] = 1 - cos
        return sintheta, costheta, multiplier, is_zero_angle

    def partial_calculate(self, group, cache):
        '''
        :param group: the group that was just permuted. repreents the indexes in self.perm that need to have A,B calculated
        '''
        # permuted_array=array[perm[i]]
        # perms[i] = [perm[perms[i - 1][j]] for j in range(size)]
        partial_calculate(self, group, cache)

    def close_cycle(self,group, cache):
        pc= self.PartialCalculation.copyconstruct(self._calc)
        self.partial_calculate(group, cache)
        return pc

    def unclose_cycle(self, calc):
        self._calc=calc

class MoleculeLegalPermuter:
    """
    This class builds a permutation atom by atom, checking with each atom whether its new position creates an illegal permutation
    (as defined by the permchecker class)
    To that end, the class uses a class called pip (perm in progress)
    The pip is created stage by stage-- each equivalency group is built atom-by-atom (into legal cycles)
    """

    def __init__(self, mol, op_order, op_type, permchecker=TruePermChecker, pipclass=CABPermInProgress):
        self._perm_count = 0
        self._groups = mol.equivalence_classes
        self._pip = pipclass(mol, op_order, op_type, permchecker)
        print(self._pip.type)
        self._cycle_lengths = (1, op_order)
        if op_type == 'SN':
            self._cycle_lengths = (1, 2, op_order)
        self._max_length = op_order
        self.cache=mol.cache

    def cycle_recursive_permute(self, pip, curr_atom, cycle_head, built_cycle, cycle_length, remainder):
            """
            Genereates the cycles recursively
            :param pip:  Permutation in Progress
            :param curr_atom: Next atom to add to the cycle
            :param cycle_head: The first (and last) atom of the cycle
            :param cycle_length: Length of cycle
            :param remainder: The free atoms in the group
            :return: Yields permutations (PermInProgresses)

            To start the recursion, current_atom and cycle_head are the same, meaning we have a cycle of length 1
            curr_atom<---curr_atom
            """
            # Check if this can be a complete cycle
            if cycle_length in self._cycle_lengths:
                # Yes it can, attempt to close it
                if pip.switch(curr_atom, cycle_head):  # complete the cycle (close ends of necklace)
                    built_cycle.append(cycle_head)
                    saved_calc=pip.close_cycle(built_cycle, self.cache)
                    if not remainder:  # perm has been completed
                        yield pip
                    else:
                        # cycle has been completed, start a new cycle with remaining atoms
                        # As explained below, the first atom of the next cycle can be chosen arbitrarily
                        yield from self.cycle_recursive_permute(pip=pip, curr_atom=remainder[0], cycle_head=remainder[0], built_cycle=list(), cycle_length=1, remainder=remainder[1:])
                    pip.unclose_cycle(saved_calc)
                    built_cycle.remove(cycle_head)
                    pip.unswitch(curr_atom, cycle_head)  # Undo the last switch
            # We now have a partial cycle of length cycle_length (we already checked it as a full cycle
            # above), now we try to extend it
            if cycle_length < self._max_length:
                for next_atom in remainder:
                    # Go over all the possibilities for the next atom in the cycle
                    if pip.switch(curr_atom, next_atom):
                        next_remainder = list(remainder)
                        next_remainder.remove(next_atom)
                        built_cycle.append(next_atom)
                        yield from self.cycle_recursive_permute(pip, next_atom, cycle_head, built_cycle, cycle_length + 1, next_remainder)
                        built_cycle.remove(next_atom)
                        pip.unswitch(curr_atom, next_atom)

    def _group_permuter(self, group, pip):
        # Start the recursion. It doesn't matter which atom is the first in the cycle, as the cycle's starting points are\
        # meaningless: 1<--2, 2<--3, 3<--1 is the same as 2<--3, 3<--1, 1<--2.
        yield from self.cycle_recursive_permute(pip=pip, curr_atom=group[0], cycle_head=group[0], built_cycle=list(), cycle_length=1, remainder=group[1:])

    def group_recursive_permute(self,groups, pip):
            if not groups:
                self._perm_count += 1
                yield pip
            else:
                for perm in self._group_permuter(groups[0], pip):
                    #saved_calc=pip.close_cycle(groups[0], self.cache)
                    yield from self.group_recursive_permute(groups[1:], perm)
                    #pip.unclose_cycle(saved_calc)

    def permute(self):
        # permutes molecule by groups
        for pip in self.group_recursive_permute(self._groups, self._pip):
            yield pip

