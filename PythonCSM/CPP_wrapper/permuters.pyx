# # cython: profile=True
# cython: language-level=3
# cython: boundscheck=False, wraparound=False, nonecheck=False
import ctypes

include "misc.pxi"
include "cache.pyx"
include "fast_calculations.pyx"
import numpy as np
cimport numpy as np
cimport cython
import math
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy


__author__ = 'Devora'


cdef class PermChecker:
    def __init__(self):
        pass

    def is_legal(self, pip, origin, destination):
        raise NotImplemented

cdef class TruePermChecker(PermChecker):
    def __init__(self, mol):
        pass

    def is_legal(self, pip, origin, destination):
        return True

cdef class LegalPermChecker(PermChecker):
    def __init__(self, mol):
        self.mol = mol

    def is_legal(self, pip, origin, destination):
        for adjacent in self.mol.atoms[destination].adjacent:
            if pip.p[adjacent] != -1 and (origin, pip.p[adjacent]) not in self.mol.bondset:
                return False
        return True

cdef class PermsHolder:
    cdef long *buffer
    cdef int molecule_size
    cdef int op_order
    cdef int buf_size

    def __cinit__(PermsHolder self, int molecule_size, int op_order):
        self.op_order = op_order
        self.molecule_size = molecule_size
        self.buf_size = op_order * molecule_size * sizeof(long)
        self.buffer = <long *>malloc(self.buf_size)

    def __dealloc__(PermsHolder self):
        if self.buffer:
            free(self.buffer)
            self.buffer = <long *>0

    cdef public PermsHolder copy(PermsHolder self):
        copy = PermsHolder(self.molecule_size, self.op_order)
        memcpy(copy.buffer, self.buffer, self.buf_size)
        return copy

    cdef inline _get_index(PermsHolder self, int op_order, int offset):
        return (op_order * self.molecule_size + offset)

    cdef public inline long get_perm_value(PermsHolder self, int op_order, int offset):
        return self.buffer[self._get_index(op_order, offset)]

    cdef public inline void set_perm_value(PermsHolder self, int op_order, int offset, int value):
        self.buffer[self._get_index(op_order, offset)] = value

    def _check_indices(PermsHolder self, indices):
        if indices[0]<0 or indices[0] >= self.op_order:
            raise ValueError("op_order index out of range")
        if indices[1]<0 or indices[1] >= self.molecule_size:
            raise ValueError("offset out of range")

    def __getitem__(PermsHolder self, indices):
        self._check_indices(indices)
        return self.get_perm_value(indices[0], indices[1])

    def __setitem__(PermsHolder self, indices, value):
        self._check_indices(indices)
        self.set_perm_value(indices[0], indices[1], value)

    def set_perm(PermsHolder self, int op_order, long[:] perm):
        if len(perm)!=self.molecule_size:
            raise ValueError("Wrong length of perm")
        if op_order < 0 or op_order >= self.op_order:
            raise ValueError("op_+order out of range")

        cdef int i
        for i in range(len(perm)):
            self.set_perm_value(op_order, i, perm[i])


    def get_perm(PermsHolder self, int op_order):
        """
        Returns a permutation as a Python array
        """
        if op_order<0 or op_order>=self.op_order:
            raise ValueError("op_order out of range")
        cdef array.array a = int_array(self.molecule_size)

        cdef int i
        for i in range(self.molecule_size):
            a[i] = self.get_perm_value(op_order, i)

        return a

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

cdef class CythonPIP:
    cdef PermChecker permchecker
    cdef CalcState _state
    cdef public str type
    cdef public int[:] p
    cdef int _size
    cdef int op_order
    cdef public double[:] costheta
    cdef public double[:] sintheta
    cdef double[:] multiplier
    def __init__(self, mol, op_order, op_type, permchecker):
        self.type="AB_Cython"
        self.permchecker=permchecker(mol)
        self._size=len(mol.atoms)
        self._state=CalcState(len(mol.atoms), op_order, True)
        self.p= -1 * np.ones((self._size,), dtype=np.int)
        self.op_order=op_order
        self.sintheta, self.costheta, self.multiplier= self._precalculate(op_type, op_order)

    property perm:
        def __get__(self):
            return self.p
    property perms:
        def __get__(self):
            return []

    property B:
        def __get__(self):
            B=np.zeros((3,),dtype=np.float64)
            for i in range(3):
                B[i]=self._state.B[i]
            return B

    property A:
        def __get__(self):
            A=np.zeros((3,3,), dtype=np.float64)
            for i in range(3):
                for j in range(3):
                    A[i][j]=self._state.A[i,j]
            return A
    property CSM:
        def __get__(self):
            return self._state.CSM


    cpdef switch(self, origin, destination):
        cdef int o=origin
        cdef int d=destination
        if self.permchecker.is_legal(self, origin, destination):
            self.p[o]=d
            return True
        return False

    cpdef unswitch(self, origin, destination):
        cdef int o=origin
        cdef int d=destination
        self.p[o]=-1

    cdef _precalculate(self, op_type, op_order):
        is_improper = op_type != 'CN'
        is_zero_angle = op_type == 'CS'
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
        return sintheta, costheta, multiplier

    cpdef close_cycle(self, group, cache):
        old_state = self._state.copy()
        self.partial_calculate(group, cache)
        return old_state

    cpdef unclose_cycle(self, old_state):
        self._state=old_state

    cdef partial_calculate(self, group, cache):
        cdef int iop
        cdef int j
        cdef int index, permuted_index
        cdef double dists
        for iop in range(1, self._state.op_order):
            dists=0.0
            for j in range(len(group)):
                index = group[j]
                permuted_index= self._state.perms.get_perm_value(iop - 1, self.p[index])
                #1: permute the iopth perm in index j:
                self._state.perms.set_perm_value(iop,index, permuted_index)
                #2: A+= self.multiplier[iop] * cache.outer_product_sum(index, permuted_index)
                self._state.A.add_mul(cache.outer_product_sum(index, permuted_index), self.multiplier[iop])
                #3: B+=self.sintheta[iop]*cache.cross(index, permuted_index)
                self._state.B.add_mul(cache.cross(index, permuted_index), self.sintheta[iop])
                #4:
                dists += cache.inner_product(index, permuted_index)
            self._state.CSM += self.costheta[iop] * dists


cdef class _CythonPermuter:
    def __init__(self, mol, op_order, op_type, permchecker=TruePermChecker, pipclass=CythonPIP):
        self._groups = mol.equivalence_classes
        self._pip = CythonPIP(mol, op_order, op_type, permchecker)
        print(self._pip.type)
        self._cycle_lengths = (1, op_order)
        if op_type == 'SN':
            self._cycle_lengths = (1, 2, op_order)
        self._max_length = op_order
        self.cache=Cache(mol)

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
            yield pip
