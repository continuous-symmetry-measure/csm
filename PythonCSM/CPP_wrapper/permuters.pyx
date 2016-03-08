import itertools

import math

import numpy as np
cimport numpy as np


__author__ = 'Devora'

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

ITYPE = np.int
ctypedef np.int_t ITYPE_t


class TruePermChecker:
    def __init__(self, mol):
        pass

    def is_legal(self, pip, origin, destination):
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
    cdef np.ndarray p
    cdef int op_order
    cdef public np.ndarray sintheta, costheta, multiplier

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

    def _precalculate(self, op_type, op_order):
        is_improper = op_type != 'CN'
        is_zero_angle = op_type == 'CS'
        sintheta = np.zeros(op_order, dtype=np.float64)
        costheta = np.zeros(op_order, dtype=np.float64)
        multiplier = np.zeros(op_order, dtype=np.float64)

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
    def __init__(self, mol, op_order, op_type, permchecker):
        super().__init__(mol, op_order, op_type, permchecker)
        self._calc=self.initialConstructor(op_order,self._size)

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
