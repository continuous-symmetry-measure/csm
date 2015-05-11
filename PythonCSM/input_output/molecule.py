from collections import namedtuple
from openbabel import OBAtom, OBElementTable

__author__ = 'zmbq'

def GetAtomicMass(symbol):
    tbl = OBElementTable()
    atomicNum = tbl.GetAtomicNum(symbol)
    atom = OBAtom()
    atom.SetAtomicNum(atomicNum)
    atom.SetIsotope(0)
    return atom.GetAtomicMass()

class Atom:
    def __init__(self, symbol, pos, useMass=True):
        self.symbol = symbol
        self.adjacent = []
        self.pos = pos
        if useMass and symbol != 'XX':
            self.mass = GetAtomicMass(symbol)
        else:
            self.mass = 0.0

    @property
    def mass(self):
        return self._mass

    @property
    def pos(self):
        return self._pos

    @property
    def symbol(self):
        return self._symbol


