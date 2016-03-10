# cython: profile=True
from copy import deepcopy

import numpy as np
cimport numpy as np
import math

class Factory:
    def matrix(self):
        return np.zeros((3,3,), dtype=np.float64)

    def vector(self):
        return np.zeros((3,), dtype=np.float64)

    def perms(self, perm_size, num_perms):
        return np.zeros((num_perms, perm_size, ), dtype=np.int64)


class Cache:
    def __init__(self, size, factory):
        self._matrices = []
        self._vectors = []
        self._counter = 0
        self._factory = factory
        self.cosines = np.zeros((size, ), dtype=np.float64)
        self.sines = np.zeros((size,), dtype=np.float64)

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


class CalcState:
    def __init__(self, molecule_size, op_order, factory):
        self.A = factory.matrix()
        self.B = factory.vector()
        self.perms = factory.perms(molecule_size, op_order)
        self.op_order = op_order
        self.molecule_size = molecule_size

    def __deepcopy__(self, memo):
        copy = CalcState.__new__(CalcState)
        copy.A = np.copy(self.A)
        copy.B = np.copy(self.B)
        copy.op_order = self.op_order
        copy.molecule_size = self.molecule_size
        copy.perms = np.copy(self.perms)

        return copy

    def __copy__(self):
        copy = CalcState.__new__(CalcState)
        copy.A = self.A
        copy.B = self.B
        copy.op_order = self.op_order
        copy.molecule_size = self.molecule_size
        copy.perms = self.perms


def one_iter(state, group, cache):
    cdef int i, j
    for i in range(1, state.op_order):
        for j in range(len(group)):
            from_index = group[j]
            to_index = state.perms[i-1, state.perms[0, from_index]]
            state.perms[i, from_index] = to_index
            state.A += cache.cosines[j] * cache.get_matrix(j)
            state.B += cache.sines[j] * cache.get_vector(j)
