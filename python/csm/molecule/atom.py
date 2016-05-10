from collections import namedtuple
from openbabel import OBAtom, OBElementTable

__author__ = 'zmbq'

_tbl = OBElementTable()

ChainedPermutation = namedtuple('ChainedPermutation', ['chain_perm', 'atom_perm'])

def GetAtomicMass(symbol):
    """ Returns the atomic mass of an element
    :param symbol: Element's symbol
    :return: The atomic mass
    """
    atomicNum = _tbl.GetAtomicNum(symbol)
    atom = OBAtom()
    atom.SetAtomicNum(atomicNum)
    atom.SetIsotope(0)
    return atom.GetAtomicMass()


def GetAtomicSymbol(atomic_num):
    """ Returns the element symbol from the atomic number
    :param atomic_num: Atomic number
    :return: Element symbol
    """
    return _tbl.GetSymbol(atomic_num)


class Atom:
    """ A single atom, alogn with its position and neighbors
    """
    def __init__(self, symbol, pos, useMass=True, chain=''):
        self._symbol = symbol
        self.adjacent = []
        self.pos = pos
        if useMass and symbol != 'XX':
            self._mass = GetAtomicMass(symbol)
        else:
            self._mass = 1.0
        self._chain = chain
        self._equivalency=[]

    @property
    def mass(self):
        return self._mass

    @property
    def symbol(self):
        return self._symbol

    @property
    def chain(self):
        return self._chain

    @property
    def equivalency(self):
        return self._equivalency

    def add_equivalence(self, index):
        self._equivalency.append(index)

    def __str__(self):
        return "Symbol: %s\tPos: %s\tAdjacent: %s" % (self.symbol, self.pos, self.adjacent)
