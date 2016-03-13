import math
import numpy as np

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

    @staticmethod
    cdef buffer_copy(double buf[3][3]):
        cdef Matrix3D copy = Matrix3D()
        cdef int i, j
        for i in range(3):
            for j in range(3):
                copy.buf[i][j] = buf[i][j]
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

    def __str__(self):
        arr=np.zeros((3,3,))
        for i in range(3):
            for j in range(3):
                arr[i][j]=self.buf[i][j]
        return str(arr)


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

    @staticmethod
    cdef buffer_copy(double buf[3]):
        cdef Vector3D copy = Vector3D()
        cdef int i
        for i in range(3):
            copy.buf[i] = buf[i]
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

    def __str__(self):
        arr=np.zeros((3,))
        for i in range(3):
                arr[i]=self.buf[i]
        return str(arr)

def cross_product(a, b):
    '''
    :param a: length 3 vector
    :param b: length 3 vector
    :return: length 3 vector, cross product of a and b
    '''
    cdef double out[3]
    out[0] = a[1] * b[2] - a[2] * b[1]
    out[1] = a[2] * b[0] - a[0] * b[2]
    out[2] = a[0] * b[1] - a[1] * b[0]
    return Vector3D.buffer_copy(out)

def inner_product(a,b):
    '''
    :param a: length 3 vector
    :param b: length 3 vector
    :return: single number, inner product of a and b
    '''
    cdef double res= a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
    return res


def outer_product_sum(a, b):
    '''
    :param a: length 3 vector
    :param b: length 3 vector
    :return: 3 x3 matrix, the outer sum of a and b plus the outer sum of b and a
    '''
    cdef double out[3][3]
    cdef int i,j
    for i in range(3):
        for j in range(3):
            out[i][j]=a[i]*b[j] + b[i]*a[j]
    return Matrix3D.buffer_copy(out)

cdef class Cache:
    cdef _cross
    cdef _outer
    cdef _inner
    def __init__(self, mol):
        size=len(mol.atoms)
        self._cross= {}
        self._outer = {}
        self._inner={}
        cdef int i, j
        for i in range(size):
            for j in range(size):
                self._cross[(i,j)]= cross_product(mol.Q[i],mol.Q[j])
                self._inner[(i,j)]= inner_product(mol.Q[i],mol.Q[j])
                self._outer[(i,j)]= outer_product_sum(mol.Q[i],mol.Q[j])

    cpdef inner_product(Cache self, int i, int j):
        return self._inner[(i,j)]

    cpdef outer_product_sum(Cache self, int i, int j):
        return self._outer[(i,j)]

    cpdef cross(Cache self, int i, int j):
        return self._cross[(i,j)]
