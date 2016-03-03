import itertools

import math

import numpy as np
from calculations.pair_cache import PairCache

__author__ = 'Devora'


class PermChecker:
    def __init__(self, mol):
        self.mol = mol

    def is_legal(self, pip, origin, destination):
        for adjacent in self.mol.atoms[destination].adjacent:
            if pip.p[adjacent] != -1 and (origin, pip.p[adjacent]) not in self.mol.bondset:
                return False
        return True


class PQPermChecker:
    def __init__(self, mol):
        self.mol = mol

    def is_legal(self, pip, origin, destination):
        for adjacent in self.mol.atoms[destination].adjacent:
            if pip.p[adjacent] != -1 and (origin, pip.p[adjacent]) not in self.mol.bondset:
                return False
            for adjacent in self.mol.atoms[origin].adjacent:
                if pip.q[adjacent] != -1 and (destination, pip.q[adjacent]) not in self.mol.bondset:
                    return False
        return True


class TruePermChecker:
    def __init__(self, mol):
        pass

    def is_legal(self, pip, origin, destination):
        return True


class PQPermInProgress:
    def __init__(self, mol, op_order, op_type, permchecker=PQPermChecker):
        size = len(mol.atoms)
        self.p = [-1] * size
        self.q = [-1] * size
        self.permchecker = permchecker(mol)

    @property
    def perm(self):
        return self.p

    def switch(self, origin, destination):
        if self.permchecker.is_legal(self, origin, destination):
            assert self.p[origin] == -1 and self.q[destination] == -1
            self.p[origin] = destination
            self.q[destination] = origin
            return True
        return False

    def unswitch(self, origin, destination):
        assert self.p[origin] == destination and self.q[destination] == origin
        self.p[origin] = -1
        self.q[destination] = -1

    def close_cycle(self):
        return 1

    def unclose_cycle(self, ABPip):
        pass

class ABPartial:
    def __init__(self, ABPip):
        self.perms=np.copy(ABPip.perms)
        self.A=np.copy(ABPip.A)
        self.B=np.copy(ABPip.B)

class ABPermInProgress:
    def __init__(self, mol, op_order, op_type, permchecker=PQPermChecker):
        size = len(mol.atoms)
        self.p = [-1] * size
        self.q = [-1] * size
        self.permchecker = permchecker(mol)
        self.perms = np.empty([op_order, size], dtype=np.int)
        self.perms[0] = [i for i in range(size)]  # identity perm
        self.A = np.zeros((3, 3,))
        self.B = np.zeros((3,1,))
        self.mol = mol
        self.cache=PairCache(mol)
        self.op_order = op_order
        self.sintheta, self.costheta, self.multiplier, self.is_zero_angle = self.precalculate(op_type, op_order)
        # these are all used to calculate intermediate values:
        self._current_cycle = list()

    @property
    def perm(self):
        return self.p

    def switch(self, origin, destination):
        if self.permchecker.is_legal(self, origin, destination):
            assert self.p[origin] == -1 and self.q[destination] == -1
            self.p[origin] = destination
            self.q[destination] = origin
            self._current_cycle.append(destination)  # check that it's this and not origin
            return True
        return False

    def unswitch(self, origin, destination):
        assert self.p[origin] == destination and self.q[destination] == origin
        self.p[origin] = -1
        self.q[destination] = -1
        self._current_cycle.pop()  # check that this is correct syntax

    def close_cycle(self):
        ABPip=ABPartial(self)
        self.partial_calc_AB()
        return ABPip

    def unclose_cycle(self, ABPip):
        self.perms=np.copy(ABPip.perms)
        self.A=np.copy(ABPip.A)
        self.B=np.copy(ABPip.B)
        pass

    def precalculate(self, op_type, op_order):
        is_improper = op_type != 'CN'
        is_zero_angle = op_type == 'CS'
        multiplicand = 2 * math.pi / op_order
        costheta = np.zeros(op_order)
        sintheta = np.zeros(op_order)
        multiplier = np.zeros(op_order)

        for i in range(1, op_order):
            x = 0
            if not is_zero_angle:
                x = multiplicand * i
            costheta[i] = math.cos(x)
            sintheta[i] = math.sin(x)
            if is_improper and (i % 2):
                multiplier[i] = -1 - costheta[i]
            else:
                multiplier[i] = 1 - costheta[i]
        return sintheta, costheta, multiplier, is_zero_angle

    def perm_cycle_cache(self):
        cycle = self._current_cycle
        # to permuate an array
        # permuted_array[index]=array[perm[index]]
        cycle_perm = [-1] * len(cycle)
        basis_perm = cycle

        for i in range(1, self.op_order):
            j = 0  # the cycle indices go from 0 to len cycle
            for index in (cycle):  # whereas the values within the cycle can be any number within size of molecule
                cycle_perm[j] = basis_perm[cycle[j]]  # hence, the mini cycle is built from 0 to cycle size
                self.perms[i][index] = basis_perm[cycle[j]]  # whereas the overall perm is built using the actual index value
                j += 1
            #yield cycle_perm  # this mini-cycle is used in calc_AB

    def partial_calc_AB(self):
        self.perm_cycle_cache()
        #opi=1 #the index of the op order
        #for cycle_perm in self.perm_cycle_cache():  # there are range(1, op_order) perms
        #    for j in range(len(cycle_perm)):
        #        k=cycle_perm[j] #the index we're dealing with
        #        self._temp_A += self.multiplier[opi] * self.cache.outer_product_sum(k, cycle_perm[i])
        #        self._temp_B += self.sintheta[opi] * self.cache.cross(k, cycle_perm[i])
        #        opi += 1
        #self.A += self._temp_A
        #self.B += self._temp_B


class MoleculeLegalPermuter:
    """
    This class builds a permutation atom by atom, checking with each atom whether its new position creates an illegal permutation
    (as defined by the permchecker class)
    To that end, the class uses a class called pip (perm in progress)
    The pip is created stage by stage-- each equivalency group is built atom-by-atom (into legal cycles)
    """

    def __init__(self, mol, op_order, op_type, permchecker=PQPermChecker, pipclass=ABPermInProgress):
        self._perm_count = 0
        self._groups = mol.equivalence_classes
        self._pip = pipclass(mol, op_order, op_type, permchecker)
        self._cycle_lengths = (1, op_order)
        if op_type == 'SN':
            self._cycle_lengths = (1, 2, op_order)
        self._max_length = op_order

    def _group_permuter(self, group, pip):
        """
        Generates permutations with cycles of a legal sizes
        """

        def recursive_permute(pip, curr_atom, cycle_head, cycle_length, remainder):
            """
            Genereates the cycles recursively
            :param pip:  Permutation in Progress
            :param curr_atom: Next atom to add to the cycle
            :param cycle_head: The first (and last) atom of the cycle
            :param cycle_length: Length of cycle
            :param remainder: The free atoms in the group
            :return: Yields permutations (PermInProgresses)

            To start the recursion, current_atom and cycle_head are the same, meaning we have a cycle of length 1
            curr_atom<---curr_atom
            """

            # Check if this can be a complete cycle
            if cycle_length in self._cycle_lengths:
                # Yes it can, close it
                if pip.switch(curr_atom, cycle_head):  # complete the cycle (close ends of necklace)
                    ABPip=pip.close_cycle()
                    if not remainder:  # perm has been completed
                        yield pip
                    else:
                        # cycle has been completed, start a new cycle with remaining atoms
                        # As explained below, the first atom of the next cycle can be chosen arbitrarily
                        yield from recursive_permute(pip, remainder[0], remainder[0], 1, remainder[1:])
                    pip.unclose_cycle(ABPip)
                    pip.unswitch(curr_atom, cycle_head)  # Undo the last switch

            # We now have a partial cycle of length cycle_length (we already checked it as a full cycle
            # above), now we try to extend it
            if cycle_length < self._max_length:
                for next_atom in remainder:
                    # Go over all the possibilities for the next atom in the cycle
                    if pip.switch(curr_atom, next_atom):
                        next_remainder = list(remainder)
                        next_remainder.remove(next_atom)
                        yield from recursive_permute(pip, next_atom, cycle_head, cycle_length + 1, next_remainder)
                        pip.unswitch(curr_atom, next_atom)

        # Start the recursion. It doesn't matter which atom is the first in the cycle, as the cycle's starting points are\
        # meaningless: 1<--2, 2<--3, 3<--1 is the same as 2<--3, 3<--1, 1<--2.
        yield from recursive_permute(pip, group[0], group[0], 1, group[1:])

    def permute(self):
        # permutes molecule by groups
        def recursive_permute(groups, pip):
            if not groups:
                self._perm_count += 1
                yield pip
            else:
                for perm in self._group_permuter(groups[0], pip):
                    yield from recursive_permute(groups[1:], perm)

        for pip in recursive_permute(self._groups, self._pip):
            yield pip


class SinglePermPermuter:
    """ A permuter that returns just one permutation, used for when the permutation is specified by the user """

    class SinglePermInProgress:
        def __init__(self, perm):
            self.perm = perm

    def __init__(self, perm):
        self._perm = self.SinglePermInProgress(perm)

    def permute(self):
        yield self._perm
