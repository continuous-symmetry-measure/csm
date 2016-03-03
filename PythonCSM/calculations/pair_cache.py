def cross_product(a, b):
    '''
    :param a: length 3 vector
    :param b: length 3 vector
    :return: length 3 vector, cross product of a and b
    '''
    out=[0,0,0]
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
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]


def outer_product_sum(a, b):
    '''
    :param a: length 3 vector
    :param b: length 3 vector
    :return: 3 x3 matrix, the outer sum of a and b plus the outer sum of b and a
    '''
    out=[[0,0,0],[0,0,0],[0,0,0]]
    for i in range(3):
        for j in range(3):
            out[i][j]=a[i]*b[j] + b[i]*a[j]
    return out


class PairFacts:
    def __init__(self,inner_product, outer_product, cross):
        self._inner=inner_product
        self._outer=outer_product
        self._cross=cross
    @property
    def cross(self):
        return self._cross

    @property
    def inner(self):
        return self._inner

    @property
    def outer(self):
        return self._outer


class PairCache:
    def __init__(self,mol):
        self.mol=mol
        self.pairfacts={}

    def calc_i_j(self, i, j):
        cross= cross_product(self.mol.Q[i],self.mol.Q[j])
        inner=inner_product(self.mol.Q[i],self.mol.Q[j])
        outer=outer_product_sum(self.mol.Q[i],self.mol.Q[j])
        self.pairfacts[(i,j)]=PairFacts(inner, outer, cross)

    def inner_product(self, i,j):
        if (i,j) in self.pairfacts:
            return self.pairfacts[(i,j)].inner
        self.calc_i_j(i,j,)
        return self.pairfacts[(i,j)].inner

    def outer_product_sum(self, i,j):
        if (i,j) in self.pairfacts:
            return self.pairfacts[(i,j)].outer
        self.calc_i_j(i,j,)
        return self.pairfacts[(i,j)].outer

    def cross(self,i,j):
        if (i,j) in self.pairfacts:
            return self.pairfacts[(i,j)].cross
        self.calc_i_j(i,j,)
        return self.pairfacts[(i,j)].cross
