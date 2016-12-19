import numpy as np
from csm.fast import PreCalcPIP

__author__ = 'Devora'

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
        ic.constraints=np.array([np.array(self.constraints[i]) for i in range(self.len)])
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
        self.choose() #debug

    def remove_constraint_from_index(self, index, constraint):
        self.constraints[index][constraint]=0
        debug=None

    def remove_index(self, index):
        try:
            self.allowed_indices.remove(index)
        except:
            pass

    def check(self):
        for index in self.allowed_indices:
            if np.sum(self.constraints[index])==0:
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
        if len(self.allowed_indices)<1:
            return None, None
        sums= [np.sum(self.constraints[i]) for i in self.allowed_indices]
        pre_index=np.argmin(sums)
        min_index=self.allowed_indices[pre_index]

        return min_index, [i for i, val in enumerate(self.constraints[min_index]) if val==1]





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

