# cython: profile=True
# cython: language-level=3
# cython: boundscheck=False, wraparound=False

include "misc.pxi"

from copy import deepcopy

import numpy as np
cimport numpy as np
cimport cython
import math

cdef class Matrix3D:
    cdef double buf[3][3]

    @staticmethod
    def zero():
        cdef Matrix3D m = Matrix3D()
        cdef int i, j
        for i in range(3):
            for j in range(3):
                m[i][j] = 0
        return m

    cdef copy(Matrix3D self):
        cdef Matrix3D copy = Matrix3D()
        cdef int i, j
        for i in range(3):
            for j in range(3):
                copy.buf[i][j] = self.buf[i][j]
        return copy

    def __getitem__(Matrix3D self, indices):
        return self.buf[indices[0]][indices[1]]

    def __setitem__(Matrix3D self, indices, v):
        self.buf[indices[0]][indices[1]] = v

    def __add__(Matrix3D first, Matrix3D second):
        cdef Matrix3D result = Matrix3D()
        cdef int i, j

        for i in range(3):
            for j in range(3):
                result.buf[i][j] = first.buf[i][j] + second.buf[i][j]
        return result

    def __mul__(Matrix3D matrix, double scalar):
        cdef Matrix3D result = Matrix3D()
        cdef int i, j
        for i in range(3):
            for j in range(3):
                result.buf[i][j] = matrix.buf[i][j] * scalar
        return result

    def __iadd__(Matrix3D self, Matrix3D other):
        cdef int i, j
        for i in range(3):
            for j in range(3):
                self.buf[i][j] += other.buf[i][j]
        return self

    def __imul__(self, scalar):
        cdef int i, j
        for i in range(3):
            for j in range(3):
                self.buf[i][j] *= scalar
        return self


cdef class Vector3D:
    cdef double buf[3]

    @staticmethod
    def zero():
        cdef Vector3D v = Vector3D()
        cdef int i
        for i in range(3):
            v[i] = 0
        return v

    cdef copy(Vector3D self):
        cdef Vector3D copy = Vector3D()
        cdef int i
        for i in range(3):
            copy.buf[i] = self.buf[i]
        return copy

    def __getitem__(Vector3D self, index):
        return self.buf[index]

    def __setitem__(Vector3D self, index, value):
        self.buf[index] = value

    def __add__(Vector3D first, Vector3D second):
        cdef Vector3D result = Vector3D()
        cdef int i

        for i in range(3):
            result.buf[i] = first.buf[i] + second.buf[i]
        return result

    def __mul__(Vector3D vector, double scalar):
        cdef Vector3D result = Vector3D()
        cdef int i
        for i in range(3):
            result.buf[i] = vector.buf[i] * scalar
        return result

    def __iadd__(Vector3D self, Vector3D other):
        cdef int i
        for i in range(3):
            self.buf[i] += other.buf[i]
        return self

    def __imul__(self, scalar):
        cdef int i
        for i in range(3):
            self.buf[i] *= scalar
        return self


class Factory:
    def matrix(self):
        # return np.zeros((3,3,), dtype=np.double)
        return Matrix3D()

    def vector(self):
        # return np.zeros((3,), dtype=np.double)
        return Vector3D()

    def perms(self, perm_size, num_perms):
        return np.zeros((num_perms, perm_size, ), dtype=np.int64)



class Cache:
    def __init__(self, size, factory):
        self._matrices = []
        self._vectors = []
        self._counter = 0
        self._factory = factory
        self.cosines = np.zeros((size, ), dtype=np.double)
        self.sines = np.zeros((size,), dtype=np.double)

        angle_inc = 2 * math.pi / float(size)
        angle = 0.0
        for i in range(size):
            self._matrices.append(self.create_matrix())
            self._vectors.append(self.create_vector())
            self.cosines[i] = math.cos(angle)
            self.sines[i] = math.sin(angle)
            angle += angle_inc

    def create_matrix(self):
        m = self._factory.matrix()
        for i in range(3):
            for j in range(3):
                m[i, j] = self._counter
                self._counter += 1.02
        return m

    def create_vector(self):
        v = self._factory.vector()
        for i in range(3):
            v[i] = self._counter
            self._counter -= 0.53
        return v

    def get_matrix(self, num):
        return self._matrices[num]

    def get_vector(self, num):
        return self._vectors[num]


cdef class CalcState:
    cdef public Matrix3D A
    cdef public Vector3D B
    cdef public perms
    cdef public int op_order
    cdef public int molecule_size

    cdef long *perm_buf

    def __init__(self, molecule_size, op_order, factory):
        if factory:
            self.A = factory.matrix()
            self.B = factory.vector()
            self.perms = factory.perms(molecule_size, op_order)
            self.set_pointers()
        self.op_order = op_order
        self.molecule_size = molecule_size
        self.perm_buf = <long *>0


    cdef set_pointers(CalcState self):
        cdef long ptr = self.perms.__array_interface__['data'][0]
        self.perm_buf = <long *>ptr

    def __deepcopy__(self, memo):
        copy = CalcState(self.molecule_size, self.op_order, None)
        copy.A = self.A.copy()
        copy.B = self.B.copy()
        copy.perms = np.copy(self.perms)
        copy.set_pointers()

        return copy

    cdef long *get_perm_value_ptr(CalcState self, int perm, int index):
        return &self.perm_buf[perm * self.molecule_size + index]

    cdef inline long get_perm_value(CalcState self, int perm, int index):
        return self.get_perm_value_ptr(perm, index)[0]

    cdef inline void set_perm_value(CalcState self, int perm, int index, int value):
        self.get_perm_value_ptr(perm, index)[0] = value


def one_iter(state, group, cache):
    cdef int i, j
    cdef int from_index, to_index
    for i in range(1, state.op_order):
        for j in range(len(group)):
            from_index = group[j]
            zero_index = state.get_perm_value(0, from_index)
            to_index = state.get_perm_value(i-1, zero_index)
            state.set_perm_value(i, from_index, to_index)
            state.A += cache.get_matrix(j) * cache.cosines[j]
            state.B += cache.get_vector(j) * cache.sines[j]
