import numpy as np
from csm.fast import calcstate_from_python


DTYPE=np.float64
ITYPE=np.long

class CalcState:
    """
    A class that stores the ongoing partial calculation results of
    the matrix A
    the vector B
    and the CSM value
    as well as the permutation
    """


    def __init__(self, molecule_size, op_order, allocate=True):
        self.op_order = op_order
        self.molecule_size = molecule_size
        if allocate:
            self.A = np.zeros((3,3), dtype=DTYPE, order='C')
            self.B = np.zeros(3, dtype=DTYPE, order='C')
            self.perms=[]
            identity_perm=np.array([i for i in range(molecule_size)], dtype=ITYPE, order='C')
            neg_perm=np.array([-1 for i in range(molecule_size)], dtype=ITYPE, order='C')
            self.perms.append(identity_perm)
            for i in range(1,op_order):
                self.perms.append(neg_perm)
            self.perms=np.array(self.perms)
            self.CSM=1.0


    def copy(self):
        copy = CalcState(self.molecule_size, self.op_order, False)
        copy.A = np.copy(self.A)
        copy.B = np.copy(self.B)
        copy.perms = np.copy(self.perms)
        copy.CSM=self.CSM
        return copy

    @property
    def cython_state(self):
        return calcstate_from_python(self)

    @property
    def perm(self):
        return self.perms[0]

    @perm.setter
    def perm(self, val):
        self.perms[0]=val



def cross_product(a, b):
    '''
    :param a: length 3 vector
    :param b: length 3 vector
    :return: length 3 vector, cross product of a and b
    '''
    out=np.zeros(3)
    out[0] = a[1] * b[2] - a[2] * b[1]
    out[1] = a[2] * b[0] - a[0] * b[2]
    out[2] = a[0] * b[1] - a[1] * b[0]
    return out

def inner_product(a,b):
    '''
    :param a: length 3 vector
    :param b: length 3 vector
    :return: single number, inner product of a and b
    '''
    res= a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
    return res


def outer_product_sum(a, b):
    '''
    :param a: length 3 vector
    :param b: length 3 vector
    :return: 3 x3 matrix, the outer sum of a and b plus the outer sum of b and a
    '''
    out=np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            out[i][j]=a[i]*b[j] + b[i]*a[j]
    return out

class Cache:
    """
    A class that stores the results of cross, outer, and inner products of two vectors.
    Specifically, it stores all the combinations within each equivalence class.
    These results can then be retrieved from dictionaries rather than recalculated each time.
    Cannot be used on larger molecules.
    """
    def __init__(self, mol):
        self._cross= {}
        self._outer = {}
        self._inner={}
        self._mol=mol
        for group in mol.equivalence_classes:
            for i in group:
                for j in group:
                    self._cross[(i,j)]= cross_product(mol.Q[i],mol.Q[j])
                    self._inner[(i,j)]= inner_product(mol.Q[i],mol.Q[j])
                    self._outer[(i,j)]= outer_product_sum(mol.Q[i],mol.Q[j])

    def inner_product(self, i, j):
        try:
            return self._inner[(i,j)]
        except KeyError:
            return inner_product(self._mol.Q[i],self._mol.Q[j])

    def outer_product_sum(self, i, j):
        try:
            return self._outer[(i,j)]
        except KeyError:
            return outer_product_sum(self._mol.Q[i],self._mol.Q[j])

    def cross(self, i, j):
        try:
            return self._cross[(i,j)]
        except KeyError:
            return cross_product(self._mol.Q[i],self._mol.Q[j])


class FakeCache(Cache):
    """
    A class that inherits from Cache, and can return the same calculations, but does not actually cache anything.
    Can be safely used on larger molecules.
    """
    def __init__(self, mol):
        self._mol=mol

    def inner_product(self, i, j):
        return inner_product(self._mol.Q[i],self._mol.Q[j])

    def outer_product_sum(self, i, j):
        return outer_product_sum(self._mol.Q[i],self._mol.Q[j])

    def cross(self, i, j):
        return cross_product(self._mol.Q[i],self._mol.Q[j])