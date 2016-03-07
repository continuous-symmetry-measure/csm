import numpy as np
def cross_product(a, b):
    '''
    :param a: length 3 vector
    :param b: length 3 vector
    :return: length 3 vector, cross product of a and b
    '''
    out=np.zeros((3,), dtype=np.float64, order="c")
    out[0] = a[1] * b[2] - a[2] * b[1]
    out[1] = a[2] * b[0] - a[0] * b[2]
    out[2] = a[0] * b[1] - a[1] * b[0]
    return np.array(out)

def inner_product(a,b):
    '''
    :param a: length 3 vector
    :param b: length 3 vector
    :return: single number, inner product of a and b
    '''
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]


def outer_product_sum(a, b):
    '''
    :param a: length 3 vector
    :param b: length 3 vector
    :return: 3 x3 matrix, the outer sum of a and b plus the outer sum of b and a
    '''
    out= np.zeros((3, 3,))
    for i in range(3):
        for j in range(3):
            out[i][j]=a[i]*b[j] + b[i]*a[j]
    return np.array(out)



class PairCache:
    def __init__(self,mol):
        self.mol=mol
        size=len(mol.atoms)
        self._cross=np.empty((size,size,3),dtype=np.float64)
        self._inner=np.empty((size,size),dtype=np.float64)
        self._outer=np.zeros((size,size,3,3,),dtype=np.float64)
        for i in range (size):
            for j in range(size):
                self.calc_i_j(i,j)
        hi=1

    def calc_i_j(self, i, j):
        self._cross[i][j]= cross_product(self.mol.Q[i],self.mol.Q[j])
        self._inner[i][j]=inner_product(self.mol.Q[i],self.mol.Q[j])
        self._outer[i][j]=outer_product_sum(self.mol.Q[i],self.mol.Q[j])

    def inner_product(self, i,j):
        return self._inner[i][j]

    def outer_product_sum(self, i,j):
        return self._outer[i][j]

    def cross(self,i,j):
        return self._cross[i][j]
