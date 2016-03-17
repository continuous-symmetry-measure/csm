import ctypes

import numpy as np
cimport numpy as np
cimport cython
import math
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy
from libcpp cimport bool

cdef class Cache
cdef class Matrix3D
cdef class Vector3D
cdef class PermsHolder

__author__ = 'Devora'

cdef class PermChecker:
    def __init__(self):
        pass

    cpdef bool is_legal(self, pip, origin, destination):
        raise NotImplemented

cdef class TruePermChecker(PermChecker):
    def __init__(self, mol):
        pass

    cpdef bool is_legal(self, pip, origin, destination):
        return True


class PQPermChecker(PermChecker):
    def __init__(self, mol):
        self.mol = mol

    def is_legal(self, CythonPIP pip, int origin, int destination):
        cdef int adjacent
        for adjacent in self.mol.atoms[destination].adjacent:
            if pip.p[adjacent] != -1 and (origin, pip.p[adjacent]) not in self.mol.bondset:
                return False
            for adjacent in self.mol.atoms[origin].adjacent:
                if pip.q[adjacent] != -1 and (destination, pip.q[adjacent]) not in self.mol.bondset:
                    return False
        return True

cdef class CalcState:
    cdef public Matrix3D A
    cdef public Vector3D B
    cdef public PermsHolder perms
    cdef public int op_order
    cdef public int molecule_size
    cdef public double CSM

    def __init__(self, int molecule_size, int op_order, allocate=True):
        self.op_order = op_order
        self.molecule_size = molecule_size
        cdef int i
        if allocate:
            self.A = Matrix3D()
            self.B = Vector3D()
            self.perms = PermsHolder(molecule_size, op_order)
            identity_perm=np.array([i for i in range(molecule_size)])
            neg_perm=np.array([-1 for i in range(molecule_size)])
            self.perms.set_perm(0, identity_perm)
            for i in range(1,op_order):
                self.perms.set_perm(i, neg_perm)
            self.CSM=1.0


    cpdef public CalcState copy(CalcState self):
        cdef CalcState copy = CalcState(self.molecule_size, self.op_order, None)
        copy.A = self.A.copy()
        copy.B = self.B.copy()
        copy.perms = self.perms.copy()
        copy.CSM=np.copy(self.CSM)
        return copy

    property perm:
        def __get__(self):
            return self.perms.get_perm(1)

cdef class CythonPIP:
    cdef PermChecker permchecker
    cdef public CalcState state
    cdef public int[:] p
    cdef public int[:] q
    cdef public int molecule_size
    cdef public int op_order
    cdef public double[:] costheta
    cdef public double[:] sintheta
    cdef double[:] multiplier

    def __init__(self, mol, op_order, op_type, permchecker):
        self.permchecker = permchecker(mol)
        self.molecule_size =len(mol.atoms)
        self.state = CalcState(len(mol.atoms), op_order, True)
        self.p = -1 * np.ones((self.molecule_size,), dtype=np.int)  # Numpy array, but created once per molecule so no worries.
        self.q = -1 * np.ones((self.molecule_size,), dtype=np.int)
        self.op_order=op_order
        self.sintheta, self.costheta, self.multiplier= self._precalculate(op_type, op_order)


    cpdef switch(CythonPIP self, int origin, int destination):
        if self.permchecker.is_legal(self, origin, destination):
            self.p[origin]=destination
            self.q[destination]=origin
            return True
        return False

    cpdef unswitch(CythonPIP self, int origin, int destination):
        self.p[origin] = -1
        self.q[destination] = -1

    cdef _precalculate(CythonPIP self, op_type, int op_order):
        cdef bool is_improper = op_type != 'CN'
        cdef bool is_zero_angle = op_type == 'CS'
        cdef double[:] sintheta = np.zeros(op_order)
        cdef double[:] costheta = np.zeros(op_order)
        cdef double[:] multiplier = np.zeros(op_order)
        cdef int i
        cdef double theta = 0.0
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

    cpdef close_cycle(self, group, Cache cache):
        return None

    cpdef unclose_cycle(self,  CalcState old_state):
        pass



cdef class PreCalcPIP(CythonPIP):
    cpdef close_cycle(self, group, Cache cache):
        old_state = self.state.copy()
        self.partial_calculate(group, cache)
        return old_state

    cpdef unclose_cycle(self,  CalcState old_state):
        self.state = old_state

    cdef partial_calculate(self, group, Cache cache):
        cdef int iop
        cdef int j
        cdef int index, permuted_index
        cdef double dists
        for iop in range(1, self.state.op_order):
            dists=0.0
            for j in range(len(group)):
                index = group[j]
                permuted_index= self.state.perms.get_perm_value(iop - 1, self.p[index])
                #1: permute the iopth perm in index j:
                self.state.perms.set_perm_value(iop,index, permuted_index)
                #2: A+= self.multiplier[iop] * cache.outer_product_sum(index, permuted_index)
                self.state.A.add_mul(cache.outer_product_sum(index, permuted_index), self.multiplier[iop])
                #3: B+=self.sintheta[iop]*cache.cross(index, permuted_index)
                self.state.B.add_mul(cache.cross(index, permuted_index), self.sintheta[iop])
                #4:
                dists += cache.inner_product(index, permuted_index)
            self.state.CSM += self.costheta[iop] * dists


cdef class CythonPermuter:
    cdef _groups
    cdef CythonPIP _pip
    cdef _cycle_lengths
    cdef int _max_length
    cdef Cache cache
    cdef public int count

    def __init__(self, mol, op_order, op_type, perm_checker, perm_class=PreCalcPIP):
        self.count=0
        self._groups = mol.equivalence_classes
        self._pip = perm_class(mol, op_order, op_type, perm_checker)
        self._cycle_lengths = (1, op_order)
        if op_type == 'SN':
            self._cycle_lengths = (1, 2, op_order)
        self._max_length = op_order
        self.cache = Cache(mol)

    def _group_recursive_permute(self,pip, curr_atom, cycle_head, built_cycle, cycle_length, remainder):
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
                    saved_state=pip.close_cycle(built_cycle, self.cache)
                    if not remainder:  # perm has been completed
                        yield pip
                    else:
                        # cycle has been completed, start a new cycle with remaining atoms
                        # As explained below, the first atom of the next cycle can be chosen arbitrarily
                        yield from self._group_recursive_permute(pip=pip, curr_atom=remainder[0], cycle_head=remainder[0], built_cycle=list(), cycle_length=1, remainder=remainder[1:])
                    pip.unclose_cycle(saved_state)
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
                        yield from self._group_recursive_permute(pip, next_atom, cycle_head, built_cycle, cycle_length + 1, next_remainder)
                        built_cycle.remove(next_atom)
                        pip.unswitch(curr_atom, next_atom)

    def _group_permuter(self, group, pip):
        """
        Generates permutations with cycles of a legal sizes
        """
        # Start the recursion. It doesn't matter which atom is the first in the cycle, as the cycle's starting points are\
        # meaningless: 1<--2, 2<--3, 3<--1 is the same as 2<--3, 3<--1, 1<--2.
        yield from self._group_recursive_permute(pip=pip, curr_atom=group[0], cycle_head=group[0], built_cycle=list(), cycle_length=1, remainder=group[1:])

    def _recursive_permute(self,groups, pip):
        if not groups:
            yield pip
        else:
            for perm in self._group_permuter(groups[0], pip):
                yield from self._recursive_permute(groups[1:], perm)

    def permute(self):
        for pip in self._recursive_permute(self._groups, self._pip):
            self.count+=1
            yield pip.state


class SinglePermPermuter:
    """ A permuter that returns just one permutation, used for when the permutation is specified by the user """

    class SinglePIP(CythonPIP):
        def __init__(self, mol, perm, op_order, op_type):
            super().__init__(mol, op_order, op_type, TruePermChecker)
            self.p=perm
            self.cache = Cache(mol)
            self.close_cycle(perm, self.cache)

    def __init__(self, perm, mol, op_order, op_type):
        self._perm = self.SinglePermInProgress(mol, perm, op_order, op_type)

    def permute(self):
        yield self._perm.state
