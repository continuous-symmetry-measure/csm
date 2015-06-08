from openbabel import OBAtom, OBElementTable

__author__ = 'zmbq'

_tbl = OBElementTable()

def GetAtomicMass(symbol):
    atomicNum = _tbl.GetAtomicNum(symbol)
    atom = OBAtom()
    atom.SetAtomicNum(atomicNum)
    atom.SetIsotope(0)
    return atom.GetAtomicMass()

def GetAtomicSymbol(atomic_num):
    return _tbl.GetSymbol(atomic_num)

class Atom:
    def __init__(self, symbol, pos, useMass=True):
        self._symbol = symbol
        self.adjacent = []
        self._pos = pos
        if useMass and symbol != 'XX':
            self._mass = GetAtomicMass(symbol)
        else:
            self._mass = 1.0

    @property
    def mass(self):
        return self._mass

    @property
    def pos(self):
        return self._pos

    @property
    def symbol(self):
        return self._symbol

    def __str__(self):
        return "Symbol: %s\tPos: %s\tAdjacent: %s" % (self.symbol, self.pos, self.adjacent)