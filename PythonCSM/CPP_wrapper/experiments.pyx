__author__ = 'zmbq'

from cpython cimport array

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
            l[j] = j+1
        yield l
        l = temp[:]


def changer_typed(size, count):
    cdef int i, j

    l = [0] * size

    for i in range(count):
        temp = l[:]
        for j in range(size):
            l[j] = j+1
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
            ar_buf[j] = j+1
        yield ar
        ar_buf[:] = temp_buf
