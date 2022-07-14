from collections import namedtuple
import warnings

from openbabel import openbabel
from openbabel.openbabel import OBAtom
from openbabel.pybel import Atom as obel_atom
_tbl=None

__author__ = 'zmbq'



ChainedPermutation = namedtuple('ChainedPermutation', ['chain_perm', 'atom_perm'])

def GetAtomicNum(symbol):
    if _tbl:
        return _tbl.GetAtomicNum(symbol)
    return openbabel.GetAtomicNum(symbol)



def GetAtomicMass(symbol):
    """ Returns the atomic mass of an element
    :param symbol: Element's symbol
    :return: The atomic mass
    """
    atomicNum = GetAtomicNum(symbol)
    atom = OBAtom()
    atom.SetAtomicNum(atomicNum)
    atom.SetIsotope(0)
    return atom.GetAtomicMass()


def GetAtomicSymbol(atomic_num):
    """ Returns the element symbol from the atomic number
    :param atomic_num: Atomic number
    :return: Element symbol
    """
    if _tbl:
        return _tbl.GetSymbol(atomic_num)
    return openbabel.GetSymbol(atomic_num)


class Atom:
    """ A single atom, along with its position and neighbors
    """

    def __init__(self, symbol, pos, index, useMass=False, chain='Simulated Chain', res_num=-1):
        """
        :param symbol: atomic symbol (ie 'C', 'H', 'N')
        :param pos: position coordinates as a list of floats [1.0, 0.0, 0.0]
        :param index: integer index in molecule list of atoms
        :param useMass: boolean, when True atomic mass is used, when False mass=1.0
        :param chain: string, name of the chain the atom belongs to
        """
        self.index = index
        self._symbol = symbol.strip()  # make sure no white space in symbol
        self.atom_name = self._symbol # for pdb files
        self.adjacent = []
        self.pos = pos
        if useMass and symbol != 'XX':
            self._mass = GetAtomicMass(symbol)
        else:
            self._mass = 1.0
        self._chain = chain
        self._equivalency = []
        self.res_num = res_num 

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
        self._chain = c

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
        a = Atom(in_dict["symbol"], in_dict["pos"], in_dict["index"], chain=in_dict["chain"])
        a._mass = in_dict["mass"]
        a._adjacent = in_dict["adjacent"]
        a._equivalency = in_dict["equivalency"]
        return a

    def __getitem__(self, index):
        if index in (0, 1, 2):
            return self.pos[index]
        raise ValueError("Invalid Index")

    def __setitem__(self, index, value):
        if index in (0, 1, 2):
            self.pos[index] = value
        raise ValueError("Invalid Index")

