import copy
import itertools

import math
import numpy as np
from csm.fast import PreCalcPIP
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



class PythonPIP:
    def __init__(self, molecule, op_order, op_type, permchecker="Irrelevant", use_cache=True):
        self.cache=Cache(molecule)
        self.mol=molecule
        self.molecule_size =len(molecule)
        self.op_order = op_order
        self.p = -1 * np.ones((self.molecule_size,), dtype=ITYPE)  # Numpy array, but created once per molecule so no worries.
        self.q = -1 * np.ones((self.molecule_size,), dtype=ITYPE)
        self._state=CalcState(self.molecule_size, op_order, True)
        self.sintheta, self.costheta, self.multiplier= self._precalculate(op_type, op_order)
        #self.permchecker = permchecker(mol)
        #self.truecount=0
        #self.falsecount=0

    @property
    def state(self):
        return self._state.cython_state

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

    def switch(self, origin, destination):
        self.p[origin]=destination
        self.q[destination]=origin

    def unswitch(self, origin, destination):
        self.p[origin]=-1
        self.q[destination]=-1

    def close_cycle(self, group):
        old_state = self._state.copy()
        self.partial_calculate(group, self.cache)
        return old_state

    def unclose_cycle(self,  old_state):
        self._state = old_state

    def partial_calculate(self, group, cache):
        for iop in range(1, self._state.op_order):
            dists=0.0
            for index in group:
                permuted_index= self._state.perms[iop - 1][self.p[index]]
                #1: permute the iopth perm in index j:
                self._state.perms[iop][index]= permuted_index
                #2: A+= self.multiplier[iop] * cache.outer_product_sum(index, permuted_index)
                self._state.A+= self.cache.outer_product_sum(index, permuted_index) * self.multiplier[iop]  #.add_mul(self.cache.outer_product_sum(index, permuted_index), self.multiplier[iop])
                #3: B+=self.sintheta[iop]*cache.cross(index, permuted_index)
                self._state.B+= self.cache.cross(index, permuted_index) * self.sintheta[iop]#.add_mul(self.cache.cross(index, permuted_index), self.sintheta[iop])
                #4:
                dists += self.cache.inner_product(index, permuted_index)
            self._state.CSM += self.costheta[iop] * dists




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





ITYPE=np.int
DTYPE=np.float64

class ConstraintsBase:
    def __init__(self, molecule, for_copy=False):
        raise NotImplementedError

    def _create_constraints(self, molecule):
        raise NotImplementedError

    def copy(self):
        '''
        :return: a copy of the Constraints. entirely by value- must not contain any copying by reference.
        '''
        raise NotImplementedError

    def set_constraint(self, index, constraints):
        '''
        sets a given index to have the constraints given and only those
        :param index:
        :param constraints: an array of indices that are permitted values for the index
        :return:
        '''
        raise NotImplementedError

    def remove_constraint_from_all(self, constraint):
        '''
        remove a specific constraint as a permitted value from every index in constraints
        :param constraint:
        :return:
        '''
        raise NotImplementedError

    def remove_constraint_from_index(self, index, constraint):
        '''
        remove a specific constraint from a specific index
        :param index:
        :param constraint:
        :return:
        '''
        raise NotImplementedError

    def remove_index(self, index):
        '''
        "removing" an index from constraints means ensuring the index is inaccessible to the function choose
        :param index:
        :return:
        '''
        raise NotImplementedError

    def check(self):
        '''
        check ascertains whether a given set of constraints represents a dead end on the permutation tree
        :return:
        '''
        raise NotImplementedError


    def __getitem__(self, item):
        '''returns atom at index, and its options as an array'''
        raise NotImplementedError

    def choose(self):
        '''chooses the "best" atom and its options as an array'''
        raise NotImplementedError


class IndexConstraints(ConstraintsBase):
    '''
    A class that uses a numpy array of size NxN (N= len(molecule)) to store constraints
    '''
    def __init__(self, molecule, for_copy=False):
        if not for_copy:
            self._create_constraints(molecule)

    def _create_constraints(self, molecule):
        self.len = len(molecule)
        self.constraints=np.zeros((len(molecule), len(molecule)), dtype='int')
        self.allowed_indices=[i for i in range(len(molecule))]
        for index, atom in enumerate(molecule.atoms):
            debug=0
            for allowed_index in atom.equivalency:
                self.constraints[index][allowed_index]=1

    def copy(self):
        ic=IndexConstraints(None, True)
        ic.constraints=np.array(self.constraints)
        ic.allowed_indices=[x for x in self.allowed_indices]
        ic.len=self.len
        return ic

    def set_constraint(self, index, constraints):
        self.constraints[index]=np.zeros(self.len)
        for allowed in constraints:
            self.constraints[index][allowed]=1

    def remove_constraint_from_all(self, constraint):
        for index in self.allowed_indices:
            self.remove_constraint_from_index(index, constraint)

    def remove_constraint_from_index(self, index, constraint):
        self.constraints[index][constraint]=0
        debug=None

    def remove_index(self, index):
        try:
            self.allowed_indices.remove(index)
        except:
            pass

    def check(self):
        if len(self.allowed_indices)==0:
            return False
        for index in self.allowed_indices:
            sum= np.sum(self.constraints[index])
            if sum<1:
                return False
        return True

    def __getitem__(self, item):
        vals=self.constraints[item]
        constraints=[]
        for i in range(len(vals)):
            if vals[i]==1:
                constraints.append(i)
        return constraints

    def choose(self):
        if len(self.allowed_indices)==0:
            return None, None
        minimum=MAXDOUBLE
        min_index=-1
        for index in self.allowed_indices:
            sum= np.sum(self.constraints[index])
            if sum<minimum:
                minimum=sum
                min_index=index
        return min_index, [x for x in self.constraints[min_index]]





class DictionaryConstraints(ConstraintsBase):
    def __init__(self, molecule, for_copy=False):
        self.constraints={}
        if not for_copy:
            self._create_constraints(molecule)


    def _create_constraints(self, molecule):
        for index, atom in enumerate(molecule.atoms):
            self.constraints[index] = list(atom.equivalency)

    def copy(self):
        dc= DictionaryConstraints(None, True)
        for index in self.constraints:
            dc.constraints[index]=list(self.constraints[index])
        return dc

    def set_constraint(self, index, constraints):
        self.constraints[index]=constraints

    def remove_constraint_from_all(self, constraint):
        for key in self.constraints:
            self.remove_constraint_from_index(key, constraint)
        pass

    def remove_constraint_from_index(self, index, constraint):
        try:
            if constraint in self.constraints[index]:
                self.constraints[index].remove(constraint)
        except KeyError:
            pass

    def remove_index(self, index):
        self.constraints.pop(index)

    def check(self):
        for key in self.constraints:
            if len(self.constraints[key])<1:
                return False
        return True

    def __getitem__(self, item):
        return self.constraints[item]

    def choose(self):
        keys = []
        lengths = []
        for key in self.constraints:
            keys.append(key)
            lengths.append(len(self.constraints[key]))

        if lengths:
            index = np.argmin(lengths)
            key = keys[index]
            return key, list(self.constraints[key])

        return None, None




class ConstraintPropagator:
    def __init__(self, molecule, op_order, op_type):
        self.molecule=molecule
        self.op_order=op_order
        self.op_type=op_type


    def propagate(self, constraints, pip, origin, destination, cycle_length, cycle_head, cycle_tail, keep_structure):
        #this is the function which handles the logic of which additional constraints get added after a placement
        #depropagating, for now, is handled by copying old constraints
        
        #the logic of constraints after a placement:

        #1: handle cycles
        self.handle_cycles(constraints, pip, cycle_length, cycle_head, cycle_tail)

        #2: any index which had destination as an option, destination is removed as an option
        constraints.remove_constraint_from_all(destination)

        #(4: keep_structure)
        if keep_structure:
            self.keep_structure(constraints, origin, destination)

        #finally: origin is removed entirely from the constraints dictionary-- it has been placed, and has nothing left
        #this is the only time and only way it is legal for an atom to have no options.
        constraints.remove_index(origin)

    def handle_cycles(self, constraints, pip, cycle_length, cycle_head, cycle_tail):
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
            constraints.remove_constraint_from_all(cycle_head)
            constraints.set_constraint(cycle_tail, [cycle_head])

        else:
            # otherwise, either cycle length is NOT op_order AND not SN:
            # in which case, it is forbidden for cycle_tail to attach to cycle head:
            if cycle_length>1:#ie 2-1
                constraints.remove_constraint_from_index(cycle_tail, cycle_head)
                if self.op_type == 'SN':
                    second_index=pip.p[cycle_head]
                    constraints.remove_constraint_from_index(second_index, cycle_head) #because it was SN, this was left unhandled while
                    # possibility of len2 cycle remained
            else: #it is length two, and is not SN, remove as normal
                if self.op_type != 'SN':
                    constraints.remove_constraint_from_index(cycle_tail, cycle_head)

    def keep_structure(self, constraints, origin, destination):
        #1. atoms in A's bondset, can only go to constraints of theirs, that are in B's bondset
        for atom_in_A_bondset in self.molecule.atoms[origin].adjacent:
            try:
                new_constraints=[]
                for constraint in constraints[atom_in_A_bondset]:
                    if (destination, constraint) in self.molecule.bondset:
                        new_constraints.append(constraint)
                        #self._remove(atom_in_A_bondset, constraint)
                constraints.set_constraint(atom_in_A_bondset, new_constraints)
            except KeyError:
                pass
            except:
                print("aha!")
        #2.






class ConstraintPermuter:
    def __init__(self, molecule, op_order, op_type, keep_structure, *args, **kwargs):
        self.molecule=molecule
        self.op_order=op_order
        self.op_type=op_type
        self.keep_structure=keep_structure
        self.count=0
        self.truecount=0
        self.falsecount=0
        self.cycle_lengths=[1,op_order]
        if op_type=='SN':
            self.cycle_lengths.append(2)
        self.constraints_prop=ConstraintPropagator(self.molecule, self.op_order, self.op_type)

    def create_cycle(self, atom, pip):
        group=[]
        index=atom
        while True:
            group.append(index)
            index=pip.p[index]
            #stop condition: reached end of cycle
            if index==atom:
                break
        return group

    def permute(self):
        #step 1: create initial empty pip and qip
        use_cache=True
        if len(self.molecule)>1000:
            print("Molecule size exceeds recommended size for using caching. (Perhaps you meant to use --approx?)")
            use_cache=False
        pip=PreCalcPIP(self.molecule, self.op_order, self.op_type, use_cache=use_cache)
        #pip = PythonPIP(self.molecule, self.op_order, self.op_type)
        #step 2: create initial set of constraints
        constraints=DictionaryConstraints(self.molecule)
        #constraints=IndexConstraints(self.molecule)
        #step 3: call recursive permute
        for pip in self._permute(pip, constraints):
            #print(pip.p, pip.state.perms[1])
            self.count+=1
            yield pip.state

    def _permute(self, pip, constraints):
        # step one: from the available atoms, choose the one with the least available options:
        atom, options= constraints.choose()

        # STOP CONDITION: if there are no atoms left, the permutation has been completed. yield permutation
        # (what if permutation is illegal?)
        if atom is None:
            yield pip

        # step two:
        # for each option (opt)
        else:
            for destination in options:
                #print((atom, destination))
                # save current constraints
                old_constraints=constraints.copy()
                # propagate changes in constraints
                cycle_head, cycle_tail, cycle_length, cycle=self.calculate_cycle(pip, atom, destination)
                self.constraints_prop.propagate(constraints, pip, atom, destination, cycle_length, cycle_head, cycle_tail, self.keep_structure)
                if self.truecount>70:
                    flag="for debugging inexer"
                if constraints.check():
                    self.truecount+=1
                    # make the change to pip
                    pip.switch(atom, destination)
                    #if completed cycle, close in pip, and save old state
                    if cycle_head==cycle_tail:
                        cycle=self.create_cycle(atom, pip)
                        old_state=pip.close_cycle(cycle)

                    #yield from recursive create on the new pip and new constraints
                    yield from self._permute(pip, constraints)

                    #if completed cycle, restore old pip state
                    if cycle_head==cycle_tail:
                        pip.unclose_cycle(old_state)

                    #undo the change to
                    pip.unswitch(atom, destination)
                else:
                    #print("DEADEND")
                    self.falsecount+=1
                constraints=old_constraints


    def calculate_cycle(self, pip, origin, destination):
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

