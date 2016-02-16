import numpy
cimport numpy

def cross(a, b):
    return numpy.array([a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]).T
