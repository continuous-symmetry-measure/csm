# cython: profile=True
# cython: language-level=3
# cython: boundscheck=False, wraparound=False, nonecheck=False
import ctypes

include "misc.pxi"

cimport numpy as np
import numpy as np
cimport cython
import math
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy

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

    cdef add_mul(Matrix3D self, Matrix3D add, double mult):
        """
        The equivalent of self += add * mul
        """
        cdef int i, j
        for i in range(3):
            for j in range(3):
                self.buf[i][j] += add.buf[i][j] * mult

    def _check_indices(self, indices):
        if indices[0] < 0 or indices[0] > 2:
            raise ValueError("First index out of range")
        if indices[1] < 0 or indices[1] > 2:
            raise ValueError("Second index out of range")

    def __getitem__(Matrix3D self, indices):
        self._check_indices(indices)
        return self.buf[indices[0]][indices[1]]

    def __setitem__(Matrix3D self, indices, v):
        self._check_indices(indices)
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

    cdef add_mul(Vector3D self, Vector3D add, double mult):
        """
        The equivalent of self += add * mul
        """
        cdef int i
        for i in range(3):
            self.buf[i] += add.buf[i] * mult

    def _check_index(self, index):
        if index<0 or index>2:
            raise ValueError("index out of range")

    def __getitem__(Vector3D self, index):
        self._check_index(index)
        return self.buf[index]

    def __setitem__(Vector3D self, index, value):
        self._check_index(index)
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

    def __dealloc(PermsHolder self):
        if self.buffer:
            free(self.buffer)
            self.buffer = <long *>0

    cdef public PermsHolder copy(PermsHolder self):
        copy = PermsHolder(self.molecule_size, self.op_order)
        memcpy(copy.buffer, self.buffer, self.buf_size)
        return copy

    cdef inline get_index(PermsHolder self, int op_order, int offset):
        return (op_order * self.molecule_size + offset)

    cdef public inline long get_perm_value(PermsHolder self, int op_order, int offset):
        return self.buffer[self.get_index(op_order, offset)]

    cdef public inline void set_perm_value(PermsHolder self, int op_order, int offset, int value):
        self.buffer[self.get_index(op_order, offset)] = value

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



cdef class Cache:
    cdef _matrices
    cdef _vectors
    cdef _counter
    cdef public cosines
    cdef public sines

    def __init__(self, size):
        self._matrices = []
        self._vectors = []
        self._counter = 0
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
        m = Matrix3D()
        for i in range(3):
            for j in range(3):
                m[i, j] = self._counter
                self._counter += 1.02
        return m

    def create_vector(self):
        v = Vector3D()
        for i in range(3):
            v[i] = self._counter
            self._counter -= 0.53
        return v

    cpdef public Matrix3D get_matrix(self, num):
        return self._matrices[num]

    cpdef public Vector3D get_vector(self, num):
        return self._vectors[num]

cdef class CalcState:
    cdef public Matrix3D A
    cdef public Vector3D B
    cdef public PermsHolder perms
    cdef public int op_order
    cdef public int molecule_size

    def __init__(self, int molecule_size, int op_order, allocate=True):
        self.op_order = op_order
        self.molecule_size = molecule_size
        if allocate:
            self.A = Matrix3D()
            self.B = Vector3D()
            self.perms = PermsHolder(molecule_size, op_order)


    cpdef public CalcState copy(CalcState self):
        cdef CalcState copy = CalcState(self.molecule_size, self.op_order, None)
        copy.A = self.A.copy()
        copy.B = self.B.copy()
        copy.perms = self.perms.copy()
        return copy


cdef inline fix_perm(CalcState state, int from_index, int op_order):
    cdef int zero_perm_index = state.perms.get_perm_value(0, from_index)
    cdef int to_index = state.perms.get_perm_value(op_order-1, zero_perm_index)
    state.perms.set_perm_value(op_order, from_index, to_index)

def one_iter(CalcState state, group, Cache cache):
    cdef int i, j

    for i in range(1, state.op_order):
        for j in range(len(group)):
            fix_perm(state, group[j], i)
            state.A.add_mul(cache.get_matrix(j), cache.cosines[j])
            state.B.add_mul(cache.get_vector(j), cache.sines[j])


cdef class ArrayHolder:
    cdef array
    cdef long *ptr

    def __init__(ArrayHolder self, allocate=True):
        self.array = np.zeros((4, 12,), dtype=np.int)
        cdef long ptr = self.array.__array_interface__['data'][0]
        self.ptr = <long *>ptr

