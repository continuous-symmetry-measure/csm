import numpy as np
from csm.fast import PreCalcPIP

__author__ = 'Devora'

ITYPE=np.int
DTYPE=np.float64

class ConstraintsBase:
    '''
    this is the class that defines the functiosn expected to be present in any implementation of a constraints data structure,
    based on how the constraints are used in both the ConstraintPermuter and the ConstraintPropagator. If the decision is
    made to change something in the interface here, changes would need to be made in those classes accordingly.

    Comments explaining each function are also located here, rather than recopying them in each implementation
    '''
    def __init__(self, molecule, for_copy=False):
        raise NotImplementedError

    def _create_constraints(self, molecule):
        '''
        a function, generally called from __init__ if for_copy is false, that creates the intitial set of constraints
        from a molecule
        :param molecule:
        :return: None
        '''
        raise NotImplementedError

    def set_constraint(self, index, constraints):
        '''
        sets a given index to have the constraints given and only those. this is called in ConstraintPropagator by two
        functions: handle_cycles and keep_structure
        :param index:
        :param constraints: an array of indices that are permitted values for the index
        :return: None
        '''
        raise NotImplementedError

    def remove_constraint_from_all(self, constraint):
        '''
        remove a specific constraint as a permitted value from every index in constraints.
        called twice by ConstraintPropagator-- once to remove the destination as an option for anyone, and once in
        handle_cycles
        :param constraint:
        :return: None
        '''
        raise NotImplementedError

    def remove_constraint_from_index(self, index, constraint):
        '''
        remove a specific constraint from a specific index. used by remove_constraint_from_all,
        also used in ConstraintPropagator handle_cycles
        :param index:
        :param constraint:
        :return: None
        '''
        raise NotImplementedError

    def remove_index(self, index):
        '''
        "removing" an index from constraints means ensuring the index is inaccessible to the function choose.
        if not correctly implemented, choose may end up picking the same empty set of constraints over and over, creating
        an infinite loop...
        it is called once by ConstraintPropagator.
        :param index:
        :return: None
        '''
        raise NotImplementedError

    def check(self):
        '''
        check ascertains whether a given set of constraints represents a dead end on the permutation tree
        when check returns true, the recursivity in the permuter advances another step inwards. hence, an incorrectly
        implemented check can create infinite loops or other unpleasantness
        :return: False if constraints have reached a dead end, True otherwise
        '''
        raise NotImplementedError


    def __getitem__(self, item):
        '''
        used by keep_structure in ConstraintPropagator, and in internal uses by the various constraint implementations
        I am fairly sure this can be a return by reference, not value-- it is used for for loops.
        :param item: index for which atom to look up constraints for
        :return: returns atom's constraints as an array
        '''
        raise NotImplementedError

    def choose(self):
        '''
        called by ConstraintPermuter
        chooses the "best" atom (generally defined as the one with the least options) and its options as an array
        if there are no atoms left, it returns None, None.
        A check against None, None is used as a stop condition in ConstraintPermuter's recursion-- it should mean a
        complete permutation has been built
        :return: (atom index, options for that atom as array) OR None, None
        '''
        raise NotImplementedError

    def mark_checkpoint(self):
        """
        A Constraints object should know how to return to a previous state. This is used by the caller
        when backtracking.

        The caller calls mark_checkpoint to mark the beginning of a recursion step.
        matk_checkpoint can be called several times.
        :return: Nothing
        """
        raise NotImplementedError

    def backtrack_checkpoint(self):
        """
        Restores the state to the last saved checkpoint
        :return: Nothing
        """
        raise NotImplementedError

class IndexConstraints(ConstraintsBase):
    '''
    A class that uses a numpy array of size NxN (N= len(molecule)), filled with zeros or ones, to store constraints.
    if constraints[i][j]==1, index i is allowed to permute to index j. similarly constraints[i][j]==0 means i cannot permute to j.
    Additionally, a python list of "allowed_indices" is maintained. Only indices in this ist are allowed to be chosen by
    the function choose.

    Currently approximately 4x slower than DictionaryConstraints on molecules tested
    '''
    def __init__(self, molecule, for_copy=False):
        if not for_copy:
            self._create_constraints(molecule)

    def _create_constraints(self, molecule):
        self.len = len(molecule)
        self.constraints=np.zeros((len(molecule), len(molecule)), dtype='int')
        self.allowed_indices=[i for i in range(len(molecule))]
        for index, atom in enumerate(molecule.atoms):
            for allowed_index in atom.equivalency:
                self.constraints[index][allowed_index]=1

    def copy(self):
        ic=IndexConstraints(None, True)
        ic.constraints=np.array([np.array(self.constraints[i]) for i in range(self.len)])
        ic.allowed_indices=[x for x in self.allowed_indices]
        ic.len=self.len
        return ic

    def set_constraint(self, index, constraints):
        self.constraints[index]=np.zeros(self.len) #erase whatever was previously in the row
        for allowed in constraints:
            self.constraints[index][allowed]=1 #copy in the new constraints

    def remove_constraint_from_all(self, constraint):
        for index in self.allowed_indices:
            self.remove_constraint_from_index(index, constraint)

    def remove_constraint_from_index(self, index, constraint):
        self.constraints[index][constraint]=0

    def remove_index(self, index):
        try:
            self.allowed_indices.remove(index)
        except:
            pass

    def check(self):
        for index in self.allowed_indices:
            if np.sum(self.constraints[index])==0: #a row of zeroes only
                return False
        return True

    def __getitem__(self, item):
        return [i for i, val in enumerate(self.constraints[item]) if val==1]

        #unreachable: more verbose code saying the same thing, also a tiny bit slower
        vals=self.constraints[item]
        constraints=[]
        for i in range(len(vals)):
            if vals[i]==1:
                constraints.append(i)
        return constraints

    def choose(self):
        if len(self.allowed_indices)<1:
            return None, None
        sums= [np.sum(self.constraints[i]) for i in self.allowed_indices] #for each allowed index, count the ones- allowed indexes
        pre_index=np.argmin(sums) #find smallest amount of ones. pre_index is position within allowed_indices, not the actual index
        min_index=self.allowed_indices[pre_index] #get the actual index

        return min_index, self[min_index]


class DictionaryConstraints(ConstraintsBase):
    def __init__(self, molecule, for_copy=False):
        self.constraints={}
        if not for_copy:
            self._create_constraints(molecule)
        self.undo = []

    def _create_constraints(self, molecule):
        for index, atom in enumerate(molecule.atoms):
            self.constraints[index] = set(atom.equivalency)

    def set_constraint(self, index, constraints):
        self.push_undo('set_constraint', (index, self.constraints[index]))
        self.constraints[index] = set(constraints)

    def remove_constraint_from_all(self, constraint):
        removed_indices = set()
        for index in self.constraints:
            try:
                self.constraints[index].remove(constraint)
                removed_indices.add(index)
            except KeyError:
                pass
        #if removed_indices:
        self.push_undo('remove_constraint_from_all', (removed_indices, constraint))

    def remove_constraint_from_index(self, index, constraint):
        try:
            self.constraints[index].remove(constraint)
            self.push_undo('remove_constraint_from_index', (index, constraint))
        except KeyError:
            pass

    def remove_index(self, index):
        old_value = self.constraints.pop(index)
        if old_value is not None:
            self.push_undo('remove_index', (index, old_value))

    def check(self):
        for key in self.constraints:
            if not self.constraints[key]:
                return False
        return True

    def __getitem__(self, item):
        return list(self.constraints[item])

    def choose(self):
        min_length = 1e40
        min_key = None

        for key in self.constraints:
            key_length = len(self.constraints[key])
            if key_length < min_length:
                min_key, min_length = key, key_length

        if min_key is not None:
            return min_key, list(self.constraints[min_key])
        return None, None

    # Checkpoints
    # -----------
    # Checkpoints are implemented by a stack of instructions that restore the constraints to its previous
    # state.
    #
    # The instruction is a tuple: (instruction_name, (index, value))
    # The instructions are:
    #     set_constraint, (index, previous_constraints)
    #     remove_constraint_from_all (removed_indices, constraint)
    #     remove_constraint_from_index (index, constraint)
    #     remove_index (index, previous_constraints)
    #
    # A checkpoint is marked with a special 'mark' instruction (None argument)

    def push_undo(self, instruction, params):
        # print("push %s, %s" % (instruction, params))
        self.undo.append((instruction, params))

    def pop_undo(self):
        instruction, params = self.undo.pop()
        # print ("pop %s, %s" % (instruction, params))
        return instruction, params

    def mark_checkpoint(self):
        self.push_undo('mark', None)


    def backtrack_checkpoint(self):
        instruction, params = self.pop_undo()
        while instruction != 'mark':
            if instruction == 'set_constraint':
                self.constraints[params[0]] = params[1]
            elif instruction == 'remove_constraint_from_all':
                constraint = params[1]
                for index in params[0]:
                    self.constraints[index].add(constraint)
            elif instruction == 'remove_constraint_from_index':
                if params[0] not in self.constraints:
                    raise ValueError("Can't find %d in constraints!" % params[0])
                constraint = self.constraints[params[0]]
                constraint.add(params[1])
                # self.constraints[params[0]].add(params[1])
            elif instruction == 'remove_index':
                self.constraints[params[0]] = params[1]
            else:
                raise ValueError("Unexpected instruction %s in undo stack", instruction)

            instruction, params = self.pop_undo()


class ConstraintPropagator:
    def __init__(self, molecule, op_order, op_type):
        self.molecule = molecule
        self.op_order = op_order
        self.op_type = op_type

    def propagate(self, constraints, pip, origin, destination, cycle_length, cycle_head, cycle_tail, keep_structure):
        # This is the function which handles the logic of which additional constraints get added after a placement
        # depropagating, for now, is handled by copying old constraints
        # the logic of constraints after a placement:

        # 1: handle cycles
        self.handle_cycles(constraints, pip, cycle_length, cycle_head, cycle_tail)

        # 2: any index which had destination as an option, destination is removed as an option
        constraints.remove_constraint_from_all(destination)

        # (4: keep_structure)
        if keep_structure:
            self.keep_structure(constraints, origin, destination)

        # finally: origin is removed entirely from the constraints dictionary-- it has been placed, and has nothing left
        # this is the only time and only way it is legal for an atom to have no options.
        constraints.remove_index(origin)

    def handle_cycles(self, constraints, pip, cycle_length, cycle_head, cycle_tail):
        # This has been moved to its own function because it wasa growing long and confusing
        # NOTE: cycles of length 1 are inherently "handled" in the next two steps.
        # hence, we are only dealing with cycles of length op_order and, possibly, 2 (if SN)
        # a: if there is an atom whose destination* is origin (X->  A->B)
        # (in other words, if origin already has a value in the reverse permutation)
        if cycle_head == cycle_tail:
            return

        # at this point, there are three options:

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
        constraints = DictionaryConstraints(self.molecule)
        #constraints=IndexConstraints(self.molecule)
        #step 3: call recursive permute
        for pip in self._permute(pip, constraints):
            self.count+=1
            yield pip.state

    def _permute(self, pip, constraints):
        # step one: from the available atoms, choose the one with the least available options:
        atom, options = constraints.choose()

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
                constraints.mark_checkpoint()
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
                constraints.backtrack_checkpoint()


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

