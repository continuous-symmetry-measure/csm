import copy
import itertools

import math
import numpy as np
from csm.calculations.cython_duplicates import CalcState, Cache

from csm.calculations.constants import MAXDOUBLE

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


class TemplatePermInProgress:
    def __init__(self, mol, op_order, op_type, permchecker):
        self._size=len(mol.atoms)
        self.p=[-1] * self._size
        self.permchecker=permchecker(mol)
        self.op_order=op_order
        self.sintheta, self.costheta, self.multiplier= self._precalculate(op_type, op_order)

    @property
    def perm(self):
        return self.p

    def switch(self, origin, destination):
        if self.permchecker.is_legal(self, origin, destination):
            self.p[origin]=destination
            return True
        return False

    def unswitch(self, origin, destination):
        self.p[origin] = -1

    def _precalculate(self, op_type, op_order):
        is_improper = op_type != 'CN'
        is_zero_angle = op_type == 'CS'
        sintheta = np.zeros(op_order)
        costheta = np.zeros(op_order)
        multiplier = np.zeros(op_order)

        for i in range(1, op_order):
            if not is_zero_angle:
                theta = 2 * math.pi * i / op_order
            cos=math.cos(theta)
            costheta[i]=cos
            sintheta[i]=math.sin(theta)
            if is_improper and (i % 2):
                multiplier[i] = -1 - cos
            else:
                multiplier[i] = 1 - cos
        return sintheta, costheta, multiplier

    def close_cycle(self, group, cache):
        return None

    def unclose_cycle(self, calc):
        pass




class PQPermInProgress(TemplatePermInProgress):
    def __init__(self, mol, op_order, op_type, permchecker):
        super().__init__(mol, op_order, op_type, permchecker)
        self.q = [-1] * self._size
        self.type="PQ"

    def switch(self, origin, destination):
        if super().switch(origin, destination):
            self.q[destination] = origin
            return True
        return False

    def unswitch(self, origin, destination):
        super().unswitch(origin, destination)
        self.q[destination] = -1



class ABPermInProgress(PQPermInProgress):
    class PartialCalculation:
        def __init__(self,A,B,CSM,perms):
            self.A=np.copy(A)
            self.B=np.copy(B)
            self.CSM=CSM
            self.perms=np.copy(perms)

        @classmethod
        def copyconstruct(cls,partcalc):
            return cls(partcalc.A,partcalc.B,partcalc.CSM,partcalc.perms)

        @classmethod
        def initialConstructor(cls, op_order, size):
            perms = -1 * np.ones([op_order, size], dtype=np.int)
            perms[0] = [i for i in range(size)]  # identity perm
            A = np.zeros((3, 3,))
            B = np.zeros((3,), dtype=np.float64, order="c")
            CSM=1.0
            return cls(A,B,CSM,perms)

    def __init__(self, mol, op_order, op_type, permchecker):
        super().__init__(mol, op_order, op_type, permchecker)
        self.type="AB"
        self._state=self.PartialCalculation.initialConstructor(op_order,self._size)

    @property
    def B(self):
        return self._state.B

    @property
    def A(self):
        return self._state.A

    @property
    def perms(self):
        return self._state.perms

    @property
    def CSM(self):
        return self._state.CSM

    def _partial_calculate(self, group, cache):
        '''
        :param group: the cycle that was just permuted. represents the indexes in self.perm that need to have A,B calculated
        '''
        # permuted_array=array[perm[i]]
        # perms[i] = [perm[perms[i - 1][j]] for j in range(size)]
        for iop in range(1, self.op_order):
            dists=0.0
            for j in range(len(group)):
                index = group[j]
                permuted_index=self._state.perms[iop - 1][self.p[index]]
                self._state.perms[iop][index] = permuted_index
                self._state.A+=self.multiplier[iop] * cache.outer_product_sum(index, permuted_index)
                self._state.B+=self.sintheta[iop]*cache.cross(index, permuted_index)
                dists += cache.inner_product(index, permuted_index)
            self._state.CSM += self.costheta[iop] * dists

    def close_cycle(self,group, cache):
        pc= self.PartialCalculation.copyconstruct(self._state)
        self._partial_calculate(group, cache)
        return pc

    def unclose_cycle(self, calc):
        self._state=calc


class MoleculeLegalPermuter:
    """
    This class builds a permutation atom by atom, checking with each atom whether its new position creates an illegal permutation
    (as defined by the permchecker class)
    To that end, the class uses a class called pip (perm in progress)
    The pip is created stage by stage-- each equivalency group is built atom-by-atom (into legal cycles)
    """

    def __init__(self, mol, op_order, op_type, permchecker=TruePermChecker, pipclass=ABPermInProgress):
        self._perm_count = 0
        self._groups = mol.equivalence_classes
        self._pip = pipclass(mol, op_order, op_type, permchecker)
        print(self._pip.type)
        self._cycle_lengths = (1, op_order)
        if op_type == 'SN':
            self._cycle_lengths = (1, 2, op_order)
        self._max_length = op_order
        self.cache=PairCache(mol)

    def _group_permuter(self, group, pip):
        """
        Generates permutations with cycles of a legal sizes
        """

        def recursive_permute(pip, curr_atom, cycle_head, built_cycle, cycle_length, remainder):
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
                # Yes it can, attempt to close it
                if pip.switch(curr_atom, cycle_head):  # complete the cycle (close ends of necklace)
                    built_cycle.append(cycle_head)
                    saved_calc=pip.close_cycle(built_cycle, self.cache)
                    if not remainder:  # perm has been completed
                        yield pip
                    else:
                        # cycle has been completed, start a new cycle with remaining atoms
                        # As explained below, the first atom of the next cycle can be chosen arbitrarily
                        yield from recursive_permute(pip=pip, curr_atom=remainder[0], cycle_head=remainder[0], built_cycle=list(), cycle_length=1, remainder=remainder[1:])
                    pip.unclose_cycle(saved_calc)
                    built_cycle.remove(cycle_head)
                    pip.unswitch(curr_atom, cycle_head)  # Undo the last switch
            # We now have a partial cycle of length cycle_length (we already checked it as a full cycle
            # above), now we try to extend it
            if cycle_length < self._max_length:
                for next_atom in remainder:
                    # Go over all the possibilities for the next atom in the cycle
                    if pip.switch(curr_atom, next_atom):
                        next_remainder = list(remainder)
                        next_remainder.remove(next_atom)
                        built_cycle.append(next_atom)
                        yield from recursive_permute(pip, next_atom, cycle_head, built_cycle, cycle_length + 1, next_remainder)
                        built_cycle.remove(next_atom)
                        pip.unswitch(curr_atom, next_atom)

        # Start the recursion. It doesn't matter which atom is the first in the cycle, as the cycle's starting points are\
        # meaningless: 1<--2, 2<--3, 3<--1 is the same as 2<--3, 3<--1, 1<--2.
        yield from recursive_permute(pip=pip, curr_atom=group[0], cycle_head=group[0], built_cycle=list(), cycle_length=1, remainder=group[1:])

    def permute(self):
        # permutes molecule by groups
        def recursive_permute(groups, pip):
            if not groups:
                self._perm_count += 1
                yield pip
            else:
                for perm in self._group_permuter(groups[0], pip):
                    #saved_calc=pip.close_cycle(groups[0], self.cache)
                    yield from recursive_permute(groups[1:], perm)
                    #pip.unclose_cycle(saved_calc)

        for pip in recursive_permute(self._groups, self._pip):
            yield pip


class SinglePermPermuter:
    """ A permuter that returns just one permutation, used for when the permutation is specified by the user """

    class SinglePermInProgress(ABPermInProgress):
        def __init__(self, mol, perm, op_order, op_type):
            super().__init__(mol, op_order, op_type, TruePermChecker)
            self.p=perm
            self.cache = PairCache(mol)
            self.close_cycle(perm, self.cache)

    def __init__(self, perm, mol, op_order, op_type):
        self._perm = self.SinglePermInProgress(mol, perm, op_order, op_type)

    def permute(self):
        yield self._perm.state





class PIP:
    def __init__(self, molecule, op_type, op_order):
        self.cache=Cache(molecule)
        self.p = [-1] * len(molecule)
        self.q = [-1] * len(molecule)
        self.mol=molecule
        self.state=CalcState(len(molecule), op_order, True)
        self.sintheta, self.costheta, self.multiplier= self._precalculate(op_type, op_order)

    @property
    def cython_state(self):
        return self.state.cython_state

    def switch(self, origin, destination):
        self.p[origin]=destination
        self.q[destination]=origin

    def unswitch(self, origin, destination):
        self.p[origin]=-1
        self.q[destination]=-1

    def close_cycle(self, atom):
        old_state=self.state.copy()
        self._calculate(atom)
        return old_state

    def unclose_cycle(self, old_state):
        self.state=old_state

    def _precalculate(self, op_type, op_order):
        is_improper = op_type != 'CN'
        is_zero_angle = op_type == 'CS'
        sintheta = np.zeros(op_order)
        costheta = np.zeros(op_order)
        multiplier = np.zeros(op_order)
        theta = 0.0
        for i in range(1, op_order):
            if not is_zero_angle:
                theta = 2 * math.pi * i / op_order
            cos=math.cos(theta)
            costheta[i]=cos
            sintheta[i]=math.sin(theta)
            if is_improper and (i % 2):
                multiplier[i] = -1 - cos
            else:
                multiplier[i] = 1 - cos
        return sintheta, costheta, multiplier

    def _calculate(self, atom):
        group=[]
        index=atom
        while True:
            group.append(index)
            index=self.p[index]
            #stop condition: reached end of cycle
            if index==atom:
                break
        for iop in range(1, self.state.op_order):
            dists=0.0
            for index in group:
                permuted_index= self.state.perms[iop - 1][self.p[index]]
                #1: permute the iopth perm in index j:
                self.state.perms[iop][index]= permuted_index
                #2: A+= self.multiplier[iop] * cache.outer_product_sum(index, permuted_index)
                self.state.A+= self.cache.outer_product_sum(index, permuted_index) * self.multiplier[iop]  #.add_mul(self.cache.outer_product_sum(index, permuted_index), self.multiplier[iop])
                #3: B+=self.sintheta[iop]*cache.cross(index, permuted_index)
                self.state.B+= self.cache.cross(index, permuted_index) * self.sintheta[iop]#.add_mul(self.cache.cross(index, permuted_index), self.sintheta[iop])
                #4:
                dists += self.cache.inner_product(index, permuted_index)
            self.state.CSM += self.costheta[iop] * dists




class ConstraintManager:
    def __init__(self, molecule, op_order, op_type, constraints=None):
        self.molecule=molecule
        self.op_order=op_order
        self.op_type=op_type
        if constraints:
            self.constraints=copy.deepcopy(constraints)
        else:
            self.create_constraints()

    def create_constraints(self):
        # this is the function which creates the initial set of constraints
        #old code that might be sufficient:
        self.constraints={}
        for index, atom in enumerate(self.molecule.atoms):
            self.constraints[index]=copy.deepcopy(atom.equivalency)

    def copy(self):
        return ConstraintManager(self.molecule, self.op_order, self.op_type, self.constraints)

    def propagate(self, pip, origin, destination, cycle_length, cycle_head, cycle_tail, keep_structure):
        #this is the function which handles the logic of which additional constraints get added after a placement
        #depropagating, for now, is handled by copying old constraints
        
        #the logic of constraints after a placement:

        #1: handle cycles
        self.handle_cycles(pip, cycle_length, cycle_head, cycle_tail)

        #2: any index which had destination as an option, destination is removed as an option
        self.remove_option(destination)

        #(4: keep_structure)
        if keep_structure:
            self.keep_structure(origin, destination)

        #finally: origin is removed entirely from the constraints dictionary-- it has been placed, and has nothing left
        #this is the only time and only way it is legal for an atom to have no options.
        self.constraints.pop(origin)


    def handle_cycles(self, pip, cycle_length, cycle_head, cycle_tail):
        #This has been moved to its own function because it wasa growing long and confusing
        # NOTE: cycles of length 1 are inherently "handled" in the next two steps.
        # hence, we are only dealing with cycles of length op_order and, possibly, 2 (if SN)
        # a: if there is an atom whose destination* is origin (X->  A->B)
        # (in other words, if origin already has a value in the reverse permutation)
        if cycle_head==cycle_tail:
            return

        #at this point, there are three options:

        # if cycle_length is op_order-1, it MUST attach to cycle_head:
        if cycle_length==self.op_order-1:
            self.remove_option(cycle_head)
            self.constraints[cycle_tail] = [cycle_head]

        else:
            # otherwise, either cycle length is NOT op_order AND not SN:
            # in which case, it is forbidden for cycle_tail to attach to cycle head:
            if cycle_length>1:#ie 2-1
                self._remove(cycle_tail, cycle_head)
                if self.op_type == 'SN':
                    second_index=pip.p[cycle_head]
                    self._remove(second_index, cycle_head) #because it was SN, this was left unhandled while
                    # possibility of len2 cycle remained
            else: #it is length two, and is not SN, remove as normal
                if self.op_type != 'SN':
                    self._remove(cycle_tail, cycle_head)

    def keep_structure(self, origin, destination):
        #if origin (A) is going to destination (B)
        #THEN:
        #1. atoms in A's bondset, can only go to constraints of theirs, that are in B's bondset
        for atom_in_A_bondset in self.molecule.atoms[origin].adjacent:
            try:
                for constraint in self.constraints[atom_in_A_bondset]:
                    if (destination, constraint) not in self.molecule.bondset:
                        self._remove(atom_in_A_bondset, constraint)
            except KeyError:
                pass




    def remove_option(self, forbidden_option):
        #remove a possible destinaiton from *all* constraints
        for key in self.constraints:
            self._remove(key, forbidden_option)

    def _remove(self, key, forbidden_option):
        try:
            if forbidden_option in self.constraints[key]:
                self.constraints[key].remove(forbidden_option)
        except KeyError:
            pass

    def check(self):
        for key in self.constraints:
            if len(self.constraints[key])<1:
                return False
        return True

class SimpleAtomChooser:
    @staticmethod
    def choose(constraints):
        # the function responsible for returning an atom and its options
        # possible alternate for function responsible for identifying dead ends

        #old code that might still be fine:
        keys = []
        lengths = []
        for key in constraints:
            keys.append(key)
            lengths.append(len(constraints[key]))

        if lengths:
            index = np.argmin(lengths)
            key = keys[index]
            return key, list(constraints[key])

        return None, None


class ConstraintPermuter:
    def __init__(self, molecule, op_order, op_type, keep_structure, *args, **kwargs):
        self.molecule=molecule
        self.op_order=op_order
        self.op_type=op_type
        self.keep_structure=keep_structure
        self.count=0
        self.truecount=0
        self.falsecount=0

        #the class/algorithm used for choosing atom. may need to be associated with the class of the constraints, but
        # for now is separate
        self.atom_chooser=SimpleAtomChooser

        self.cycle_lengths=[1,op_order]
        if op_type=='SN':
            self.cycle_lengths.append(2)

    def permute(self):
        #step 1: create initial empty pip and qip
        pip=PIP(self.molecule, self.op_type, self.op_order)
        #step 2: create initial set of constraints
        constraints=ConstraintManager(self.molecule, self.op_order, self.op_type)
        #step 3: call recursive permute
        for pip in self._permute(pip, constraints):
            #print(pip.p, pip.state.perms[1])
            self.count+=1
            yield pip.cython_state

    def _permute(self, pip, constraints):
        # step one: from the available atoms, choose the one with the least available options:
        atom, options= self.atom_chooser.choose(constraints.constraints)

        # STOP CONDITION: if there are no atoms left, the permutation has been completed. yield permutation
        # (what if permutation is illegal?)
        if atom is None:
            yield pip

        # step two:
        # for each option (opt)
        else:
            for destination in options:
                # save current constraints
                old_constraints=constraints.copy()
                # propagate changes in constraints
                cycle_head, cycle_tail, cycle_length, cycle=self.calculate_cycle(pip, atom, destination, constraints)
                constraints.propagate(pip, atom, destination, cycle_length, cycle_head, cycle_tail, self.keep_structure)
                if constraints.check():
                    self.truecount+=1
                    # make the change to pip
                    pip.switch(atom, destination)
                    #if completed cycle, close in pip, and save old state
                    if cycle_head==cycle_tail:
                        old_state=pip.close_cycle(atom)

                    #yield from recursive create on the new pip and new constraints
                    yield from self._permute(pip, constraints)

                    #if completed cycle, restore old pip state
                    if cycle_head==cycle_tail:
                        pip.unclose_cycle(old_state)

                    #undo the change to
                    pip.unswitch(atom, destination)
                else:
                    self.falsecount+=1
                constraints=old_constraints


    def calculate_cycle(self, pip, origin, destination, constraints):
        cycle=[origin]
        cycle_length = 1
        index = origin
        while pip.q[index] not in [-1, origin]:
            cycle.append(pip.q[index])
            cycle_length += 1
            index = pip.q[index]

        cycle_head = index

        if cycle_head==destination:
            return cycle_head, destination, cycle_length, cycle

        # b: if destination has its own destination* (A->B ->X)
        index = destination
        while pip.p[index] not in [-1, destination, origin]:
            cycle.append(pip.p[index])
            cycle_length += 1
            index = pip.p[index]
        cycle_tail = index
        return cycle_head, cycle_tail, cycle_length, cycle





if __name__=="__main__":
    import sys
    from csm.input_output.arguments import get_split_arguments
    from csm.input_output.readers import read_inputs
    args= sys.argv[1:]
    in_args, calc_args, out_args = get_split_arguments(args)
    calc_args['molecule'], calc_args['perm'], calc_args['dirs'] = read_inputs(**in_args)
    c=ConstraintPermuter(**calc_args)

