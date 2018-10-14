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
    def __init__(self, symbol, pos, index, useMass=False, chain='Simulated Chain'):
        """
        :param symbol: atomic symbol (ie 'C', 'H', 'N')
        :param pos: position coordinates as a list of floats [1.0, 0.0, 0.0]
        :param index: integer index in molecule list of atoms
        :param useMass: boolean, when True atomic mass is used, when False mass=1.0
        :param chain: string, name of the chain the atom belongs to
        """
        self.index=index
        self._symbol=symbol.strip() #make sure no white space in symbol
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

    @chain.setter
    def chain(self, c):
        self._chain=c

    @property
    def equivalency(self):
        return self._equivalency

    def add_equivalence(self, index):
        self._equivalency.append(index)

    def __str__(self):
        return "Symbol: %s\tPos: %s\tChain: %s\tAdjacent: %s" % (self.symbol, self.pos, self.chain, self.adjacent)

    def to_dict(self):
        return {
        "index": self.index,
        "symbol": self._symbol,
        "adjacent": self.adjacent,
        "pos": self.pos,
        "mass": self._mass,
        "chain": self._chain,
        "equivalency": self._equivalency
        }

    @staticmethod
    def from_dict(in_dict):
        a=Atom(in_dict["symbol"], in_dict["pos"], in_dict["index"], chain=in_dict["chain"])
        a._mass=in_dict["mass"]
        a._adjacent=in_dict["adjacent"]
        a._equivalency=in_dict["equivalency"]
        return a


