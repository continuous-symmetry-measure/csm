__author__ = 'zmbq'

from cpython cimport array
import numpy as np
cimport numpy as np

def yielder_straight(size, count):
    l = [0] * size

    for i in range(count):
        for j in range(size):
            l[j] = j+i
        yield l

def yielder_typed(size, count):
    cdef int i, j
    l = [0] * size

    for i in range(count):
        for j in range(size):
            l[j] = j+i
        yield l


cdef array.array _int_array_template = array.array('i', [])
def yielder_array(size, count):
    cdef int i, j

    cdef array.array ar = array.clone(_int_array_template, size, False)
    cdef int[:] buf = ar

    for i in range(count):
        for j in range(size):
            buf[j] = j+i
        yield ar


def changer_straight(size, count):
    l = [0] * size

    for i in range(count):
        temp = l[:]
        for j in range(size):
            l[j] = j+i
        yield l
        l = temp[:]


def changer_typed(size, count):
    cdef int i, j

    l = [0] * size

    for i in range(count):
        temp = l[:]
        for j in range(size):
            l[j] = j+i
        yield l
        l = temp[:]


def changer_array(size, count):
    cdef int i, j

    cdef array.array ar = array.clone(_int_array_template, size, False)
    cdef array.array temp = array.clone(_int_array_template, size, False)
    cdef int[:] ar_buf = ar
    cdef int[:] temp_buf = temp

    for i in range(count):
        temp_buf[:] = ar_buf
        for j in range(size):
            ar_buf[j] = j+i
        yield ar
        ar_buf[:] = temp_buf

def second_yielder_straight(size, count):
    l = [0] * size * 2

    for half in yielder_straight(size, count):
        for i in range(size):
            l[i] = half[i] + 1
        for i in range(size):
            l[i+size] = half[i] + 1
        yield l

def second_yielder_typed(int size, int count):
    cdef int i

    l = [0] * size * 2

    for half in yielder_typed(size, count):
        for i in range(size):
            l[i] = half[i] + 1
        for i in range(size):
            l[i+size] = half[i] + 1
        yield l

def second_yielder_array(size, count):
    cdef int i
    cdef array.array ar = array.clone(_int_array_template, size*2, False)
    cdef int[:] ar_buf = ar
    cdef int[:] half_buf

    for half in yielder_array(size, count):
        half_buf = half
        for i in range(size):
            ar_buf[i] = half_buf[i] + 1
        for i in range(size):
            ar_buf[i+size] = half_buf[i] + 1
        yield ar

def measure_straight(size, count):
    container = []
    for i in range(count):
        l = [0] * size
        container.append(l)

    for l in container:
        for j in range(size):
            l[j] = j+i
    yield 0

def measure_typed(int size, int count):
    cdef int i, j
    container = []

    for i in range(count):
        l = [0] * size
        container.append(l)

    for l in container:
        for j in range(size):
            l[j] = j+i
    yield 0

def measure_array(int size, int count):
    cdef int i, j
    cdef array.array l
    container = []

    for i in range(count):
        l = array.clone(_int_array_template, size, False)
        container.append(l)

    for l in container:
        for j in range(size):
            l[j] = j+i
    yield 0

def measure_memory_view(int size, int count):
    cdef int i, j
    cdef array.array l
    cdef int[:] view

    container = []

    for i in range(count):
        l = array.clone(_int_array_template, size, False)
        container.append(l)

    for view in container:
        for j in range(size):
            view[j] = j+i
    yield 0

def measure_pointer(int size, int count):
    cdef int i, j
    cdef array.array l
    cdef int *ptr
    container = []

    for i in range(count):
        l = array.clone(_int_array_template, size, False)
        container.append(l)

    for l in container:
        ptr = l.data.as_ints
        for j in range(size):
            ptr[j] = j+i
    yield 0

DTYPE = np.int
ctypedef np.int_t DTYPE_t

cdef c_measure_numpy(int size, int count):
    cdef int i, j
    cdef np.ndarray[DTYPE_t, ndim=1] l

    container = []
    for i in range(count):
        l = np.zeros(size, dtype='int')
        container.append(l)

    for l in container:
        for j in range(size):
            l[j] = j+i

def measure_numpy(int size, int count):
    c_measure_numpy(size, count)
    yield 0
