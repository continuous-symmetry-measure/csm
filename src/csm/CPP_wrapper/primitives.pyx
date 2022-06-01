from cpython cimport array
cimport numpy as np
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy
from libcpp.vector cimport vector

np.set_printoptions(precision=20)

cdef class Matrix3D:
    cdef double buf[3][3]

    @staticmethod
    def zero():
        cdef Matrix3D m = Matrix3D()
        cdef int i, j
        for i in range(3):
            for j in range(3):
                m.buf[i][j] = 0
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
        #for i in range(3):
        #    for j in range(3):
        #        print(self.buf[i][j])
        return str(self.to_numpy())

    cpdef to_numpy(Matrix3D self):
        cdef np.ndarray[np.float_t, ndim=2] m = np.empty((3,3,), dtype=float)
        cdef int i, j
        for i in range(3):
            for j in range(3):
                m[i][j] = self.buf[i][j]
        return m

    cpdef T_mul_by_vec(Matrix3D self, Vector3D vec):
        """
        Calculates self.T @ vec
        """
        cdef Vector3D result = Vector3D.zero()
        cdef int i, j
        for i in range(3):
            for j in range(3):
                result.buf[i] += self.buf[i][j] * vec.buf[j]
        return result


cdef class Vector3D:
    cdef double buf[3]

    @staticmethod
    def zero():
        cdef Vector3D v = Vector3D()
        cdef int i
        for i in range(3):
            v.buf[i] = 0
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

    cpdef to_numpy(Vector3D self):
        cdef np.ndarray[np.float_t, ndim=1] v = np.empty((3,), dtype=float)
        cdef int i
        for i in range(3):
            v[i] = self.buf[i]
        return v


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
        #for i in range(3):
        #    print(self.buf[i])
        return str(self.to_numpy())

cdef class Vector3DHolder:
    cdef double *buffer
    cdef int buf_size

    def __cinit__(Vector3DHolder self, np.ndarray[DTYPE_t, ndim=2] mat):
        if mat.shape[1] != 3:
            raise ValueError("Expected a matrix with 3D vectors, received shape: (%d, %d)" % (mat.shape[0], mat.shape[1]))
        self.buf_size = 3 * mat.shape[0] * sizeof(double)
        self.buffer = <double *>malloc(self.buf_size)

        cdef int offset = 0
        cdef int i, j
        for i in range(mat.shape[0]):
            for j in range(3):
                self.buffer[offset] = mat[i][j]
                offset += 1

    def __dealloc__(self):
        if self.buffer:
            free(self.buffer)
            self.buffer = <double *>0

    cdef double *get_vector(self, int index):
        return &self.buffer[3 * index]



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

    property perm:
        def __get__(self):
            '''
            NOTE: this function is called in exact_calculations.py, every time we replace the best permutation with something new
            #TODO: it may be more efficient to somehow *not* reallocate this every time
            :return:
            '''
            return self.get_perm(1)

        def __set__(self, val):
                self.set_perm(1, val)



