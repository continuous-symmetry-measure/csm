# cython: profile=True

import itertools

import math

import numpy as np
cimport numpy as np

from molecule.molecule import Molecule

__author__ = 'Devora'

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

ITYPE = np.int
ctypedef np.int_t ITYPE_t


cdef class PermChecker:
    def __init__(self):
        pass

    def is_legal(self, pip, origin, destination):
        return False


class TruePermChecker(PermChecker):
    def __init__(self, mol):
        pass

    def is_legal(self, pip, origin, destination):
        return True


class LegalPermChecker(PermChecker):
    def __init__(self, mol):
        self.mol = mol

    def is_legal(self, pip, origin, destination):
        for adjacent in self.mol.atoms[destination].adjacent:
            if pip.p[adjacent] != -1 and (origin, pip.p[adjacent]) not in self.mol.bondset:
                return False
        return True


cdef class PartialCalculation:
    cdef public np.ndarray A
    cdef public np.ndarray B
    cdef public np.ndarray perms
    cdef public float CSM
    def __init__(self,np.ndarray[DTYPE_t, ndim=2]A, np.ndarray[DTYPE_t, ndim=1] B, float CSM, np.ndarray[ITYPE_t, ndim=2]perms):
        self.A=np.copy(A)
        self.B=np.copy(B)
        self.perms=np.copy(perms)
        self.CSM=np.copy(CSM)


cdef class TemplatePermInProgress:
    cdef int _size
    cdef public np.ndarray p
    cdef int op_order
    cdef public np.ndarray sintheta, costheta, multiplier
    cdef PermChecker permchecker

    def __init__(self, mol, op_order, op_type, permchecker):
        self._size=len(mol.atoms)
        self.p=np.ones(self._size, dtype=np.float64) * -1
        self.permchecker=permchecker(mol)
        self.op_order=op_order
        self.sintheta, self.costheta, self.multiplier = self._precalculate(op_type, op_order)

    property perm:
        def __get__(self):
            return self.p

    cpdef switch(self, origin, destination):
        if self.permchecker.is_legal(self, origin, destination):
            self.p[origin]=destination
            return True
        return False

    cpdef unswitch(self, origin, destination):
        self.p[origin] = -1

    cdef _precalculate(self, op_type, op_order):
        is_improper = op_type != 'CN'
        is_zero_angle = op_type == 'CS'
        sintheta = np.zeros(op_order, dtype=np.float64)
        costheta = np.zeros(op_order, dtype=np.float64)
        multiplier = np.zeros(op_order, dtype=np.float64)

        cdef int i
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
        return sintheta, costheta, multiplier

    def close_cycle(self, group, cache):
        return None

    def unclose_cycle(self, calc):
        pass


cdef class _CythonABPermInProgress(TemplatePermInProgress):
    cdef PartialCalculation _calc
    cdef public str type
    def __init__(self, mol, op_order, op_type, permchecker):
        super().__init__(mol, op_order, op_type, permchecker)
        self._calc=self.initialConstructor(op_order,self._size)
        self.type="AB_Cython"

    property CSM:
            def __get__(self):
                return self._calc.CSM
    property A:
            def __get__(self):
                return self._calc.A
    property B:
            def __get__(self):
                return self._calc.B
    property perms:
            def __get__(self):
                return self._calc.perms

    cdef copyconstruct(self,partcalc):
        return PartialCalculation(partcalc.A,partcalc.B,partcalc.CSM,partcalc.perms)

    cdef initialConstructor(self,op_order, size):
        perms = -1 * np.ones([op_order, size], dtype=np.int)
        perms[0] = [i for i in range(size)]  # identity perm
        A = np.zeros((3, 3,), dtype=np.float64)
        B = np.zeros((3,), dtype=np.float64, order="c")
        CSM=1.0
        return PartialCalculation(A,B,CSM,perms)

    cdef _partial_calculate(self, group, cache):
        '''
        :param group: the cycle that was just permuted. represents the indexes in self.perm that need to have A,B calculated
        '''
        # permuted_array=array[perm[i]]
        # perms[i] = [perm[perms[i - 1][j]] for j in range(size)]
        cdef int iop
        cdef int j
        cdef int index, permuted_index
        cdef float dists
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

    cpdef close_cycle(self,group, cache):
        pc= self.copyconstruct(self._calc)
        self._partial_calculate(group, cache)
        return pc

    cpdef unclose_cycle(self, calc):
        self._calc=calc


class CythonABPermInProgress(_CythonABPermInProgress):
    def __init__(self, mol, op_order, op_type, permchecker):
        super().__init__(mol, op_order, op_type, permchecker)
        self.type="AB_Cython"


cdef class _CythonPermuter:
    """
    This class builds a permutation atom by atom, checking with each atom whether its new position creates an illegal permutation
    (as defined by the permchecker class)
    To that end, the class uses a class called pip (perm in progress)
    The pip is created stage by stage-- each equivalency group is built atom-by-atom (into legal cycles)
    """

    def __init__(self, mol, op_order, op_type, permchecker=TruePermChecker, pipclass=None):
        self._groups = mol.equivalence_classes
        self._pip = _CythonABPermInProgress(mol, op_order, op_type, permchecker)
        print(self._pip.type)
        self._cycle_lengths = (1, op_order)
        if op_type == 'SN':
            self._cycle_lengths = (1, 2, op_order)
        self._max_length = op_order
        self.cache=mol.cache

    def _cycle_recursive_permute(self,pip, curr_atom, cycle_head, built_cycle, cycle_length, remainder):
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
                        yield from self._cycle_recursive_permute(pip=pip, curr_atom=remainder[0], cycle_head=remainder[0], built_cycle=list(), cycle_length=1, remainder=remainder[1:])
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
                        yield from self._cycle_recursive_permute(pip, next_atom, cycle_head, built_cycle, cycle_length + 1, next_remainder)
                        built_cycle.remove(next_atom)
                        pip.unswitch(curr_atom, next_atom)

    def _group_permuter(self, group, pip):
        """
        Generates permutations with cycles of a legal sizes
        """
        # Start the recursion. It doesn't matter which atom is the first in the cycle, as the cycle's starting points are\
        # meaningless: 1<--2, 2<--3, 3<--1 is the same as 2<--3, 3<--1, 1<--2.
        yield from self._cycle_recursive_permute(pip=pip, curr_atom=group[0], cycle_head=group[0], built_cycle=list(), cycle_length=1, remainder=group[1:])


    def _recursive_permute(self,groups, pip):
            if not groups:
                yield pip
            else:
                for perm in self._group_permuter(groups[0], pip):
                    #saved_calc=pip.close_cycle(groups[0], self.cache)
                    yield from self._recursive_permute(groups[1:], perm)
                    #pip.unclose_cycle(saved_calc)

    def permute(self):
        # permutes molecule by groups
        for pip in self._recursive_permute(self._groups, self._pip):
            yield pip