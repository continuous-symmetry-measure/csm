from collections import namedtuple
from openbabel import OBAtom

__author__ = 'zmbq'

_elements = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
             "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
             "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
             "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
             "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
             "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
             "Pn", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
             "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
             "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
             "Pa", "U", "Np", "Pu", "An", "Cn", "Bk", "Cf", "Es", "Fm",
             "Md", "No", "Lr", "Rf", "Ha"]

def _fillAtomicMasses():
    masses = {}
    atomicNumber = 1
    for element in _elements:
        atom = OBAtom()
        atom.SetAtomicNum(atomicNumber)
        atom.SetIsotope(0)
        masses[element] = atom.GetAtomicMass()

AtomicMasses = _fillAtomicMasses()

class Atom:
    def __init__(self, symbol, pos):
        self._symbol = symbol
        self.adjacent = []
        self._pos = pos
        self._mass = AtomicMasses[symbol]

    @property
    def mass(self):
        return self._mass

    @property
    def pos(self):
        return self._pos

    @property
    def symbol(self):
        return self._symbol

