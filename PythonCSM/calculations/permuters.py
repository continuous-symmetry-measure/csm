import itertools

__author__ = 'zmbq'


class MoleculePermuter:
    def __init__(self, mol, opOrder, is_SN):
        self._mol = mol
        self._op_order = opOrder
        self._add_cycles_of_two = is_SN
        self._circle_cache = {}  # perm_size->all circles of size
        self._CACHE_LIMIT = 10

    def _get_cycle_structs(self, perm_size, cycle_sizes):
        """
        Generates a list of cycles in a permutation. The cycles cover the entire permutation,
        and are only of sizes in cycle_sizes
        :param perm_size: Permutation size
        :param cycle_sizes: Allowed cycle sizes
        :return: A generator of a list of cycles covering the entire permtutation

        """

        def generate(cycles, left):
            len_left = len(left)
            if len_left == 0:
                yield cycles
                return

            for cycle_size in cycle_sizes:
                if cycle_size > len_left:
                    continue
                rest = left[1:]

                if cycle_size == 1:
                    yield from generate(cycles, rest)
                    continue

                start = (left[0],)
                for rest in itertools.combinations(rest, cycle_size - 1):
                    cycle = start + rest
                    new_left = [item for item in left if item not in cycle]
                    yield from generate(cycles + [cycle], new_left)

        if 1 not in cycle_sizes:
            raise ValueError("1 must be a valid cycle size")
        yield from generate([], list(range(perm_size)))

    def _calc_all_circle_permutations(self, size):
        """ Returns the permutation of the cycle """
        # To compute a full cycle of length n, we take a permutation p of size n-1 and
        # create the cycle like so: 0 goes to p[0], p[0] to p[1], p[1] to p[2] and so on, until p[n-1] goes back to 0
        # For this to work, p needs to be a permutation of 1..n-1.
        #
        # For example, size = 3 we have p=(1,2) for the circle 0->1->2->0 and p=(2,1) for the circle 0->2->1->0
        trivial = tuple(range(1, size))  # (1,2,...size-1)
        for necklace in itertools.permutations(trivial):
            # The actual necklace is [0]+necklace
            cycle_perm = [0] * size
            cur = 0
            for element in necklace:
                cycle_perm[cur] = element
                cur = element
            cycle_perm[necklace[-1]] = 0  # Add the [0] that is missing from the necklace
            yield cycle_perm

    def _all_circle_permutations(self, size):
        if size > self._CACHE_LIMIT:
            return self._calc_all_circle_permutations(size)

        if not size in self._circle_cache:
            entries = list(self._calc_all_circle_permutations(size))
            self._circle_cache[size] = entries
        result = self._circle_cache[size]
        return result

    def _all_perms_from_cycle_struct(self, perm_size, cycle_struct):
        """
        Generates all permutations with a specific cycle structure
        :param perm_size: Permutation size
        :param cycle_struct: Cycle structure (list of cycles, each cycle is itself a list of its indices)
        :return: Generator that generates all the permutations
        """
        trivial = list(range(perm_size))

        def generate(perm, cycle_index):
            # Goes over all the circles of the first cycle, apply each circle and
            # recursively generates the circles of the rest of the cycles
            # perm is built gradually, each recursive call applies one cycle
            if cycle_index < 0:
                yield tuple(perm)
                return

            cycle = cycle_struct[cycle_index]
            if len(cycle) == 1:
                yield from generate(perm, cycle_index - 1)
            else:
                # Example:
                # Lets say the permutation is (0, 1 ,2 ,3), and the cycle is (0, 1, 3)
                circles = self._all_circle_permutations(len(cycle))
                for circle_perm in circles:
                    # _all_circle_permtuations yields (1, 2, 0) and (2, 0 ,1)
                    # The permutations we need to return are (1, 3, 2, 0) and (3, 0, 2, 1) - these have
                    # one stationary point - 2, and a cycle of length 3.
                    # To do this, we need to convert the cycle from (1,2,0) and (2,0,1) to (1,3,0) and (3,0,1)
                    for i in range(len(cycle)):
                        converted_circle = cycle[circle_perm[i]]
                        perm[cycle[
                            i]] = converted_circle  # Perm is now (1, 3, 2, 0) or (3, 0, 2, 1) <--- The converted_circle applied to (0, 1, 3)
                    yield from generate(perm, cycle_index - 1)  # Apply the rest of the circles

        yield from generate(trivial[:], len(cycle_struct) - 1)

    def _group_permuter(self, group_size):
        """
        Generates all permutations of size and cycle sizes.
        :param group_size: Size of group
        :param cycle_size: Cycle size
        :param add_cycles_of_two: When true, cycles of size of two are legal
        :return: Generator for all the permutations
        """
        cycle_lengths = {1, self._op_order}
        if self._add_cycles_of_two:
            cycle_lengths.add(2)
        cycle_lengths = sorted(cycle_lengths)

        for cycle_struct in self._get_cycle_structs(group_size, cycle_lengths):
            # Loop over all cycle structures, which are of the allowed lengths covering the entire permutation
            yield from self._all_perms_from_cycle_struct(group_size,
                                                         cycle_struct)  # Return all permutations for each cycle structure

    def permute(self):
        """
        Generates all permutations of a molecule
        :param molecule_size: Molecule size
        :param groups: Equivalency groups
        :param cycle_size: Allowed cycle size
        :param add_cycles_of_two: When true, cycles of size of two are legal
        :return: Generator for all the permutations
        """

        def generate(ordered_elements, groups_left):
            """ Goes over all the permutations of the first group, applies each
             permutation and recursively aplies permutations of the entire groups
            :param ordered_elements:
            :param groups_left:
            :return:
            """
            if not groups_left:
                yield tuple(ordered_elements)
                return

            group = groups_left[0]
            groups_left = groups_left[1:]
            if len(group) == 1:
                yield from generate(ordered_elements, groups_left)
            else:
                # Example:
                # Lets say the permutation is (0, 1 ,2 ,3), and the group is (0, 1, 3)
                start_elements_order = ordered_elements[:]
                for group_perm in self._group_permuter(len(group)):
                    # group_permuter yields (1, 2, 0) and (2, 0 ,1)
                    # The permutations we need to return are (1, 3, 2, 0) and (3, 0, 2, 1) - these have
                    # one stationary point - 2, and a cycle of length 3.
                    # To do this, we need to convert the group_perm from (1,2,0) and (2,0,1) to (1,3,0) and (3,0,1)
                    converted_group = [group[i] for i in group_perm]  # Converted circle is now (1,3,0) or (3,0,1)
                    for i in range(len(converted_group)):
                        ordered_elements[group[i]] = start_elements_order[converted_group[i]]
                        # Perm is now (1, 3, 2, 0) or (3, 0, 2, 1) <--- The converted_circle applied to (0, 1, 3)
                    yield from generate(ordered_elements, groups_left)  # Apply the rest of the circles
                    ordered_elements = start_elements_order[:]

        groups = self._mol.equivalence_classes
        elements_order = [i for i in range(len(self._mol.atoms))]  # The starting elements order
        yield from generate(elements_order, groups)


class OldMoleculeLegalPermuter:
    '''
    This class builds a permutation atom by atom, checking with each atom whether its new position creates an illegal permutation (one where an atom is not connected ot its neighbors)
    To that end, the class uses a list called pip (perm in progress)
    The pip is created stage by stage-- each equivalency group is divided into cycles which are then permuted aotm by atom
    '''
    def __init__(self, mol, opOrder, is_SN):
        self._num = 0
        self._mol = mol
        self._op_order = opOrder
        self._add_cycles_of_two = is_SN
        self._empty_perm = [-1] * (len(mol.atoms))

    def _get_cycle_structs(self, group):
        """
        Generates a list of cycles in a permutation. The cycles cover the entire permutation,
        and are only of sizes in cycle_sizes
        :param perm_size: Permutation size
        :param cycle_sizes: Allowed cycle sizes
        :return: A generator of a list of cycles covering the entire permtutation

        """
        def generate(cycles, left):
            len_left = len(left)
            if len_left == 0:
                yield cycles

            for cycle_size in cycle_sizes:
                if cycle_size > len_left:
                    continue
                rest = left[1:]

                # if cycle_size==1:
                #    yield group
                #    continue

                start = (left[0],)
                for rest in itertools.combinations(rest, cycle_size - 1):
                    cycle = start + rest
                    new_left = [item for item in left if item not in cycle]
                    yield from generate(cycles + [cycle], new_left)


        cycle_sizes = {1, self._op_order}
        if self._add_cycles_of_two:
            cycle_sizes.add(2)
        return generate([], list(group))

    def _is_legal(self, pip, toi, fromi):
        # fromi,j->toi,p(j)
        for adjacent in self._mol.atoms[fromi].adjacent:
            if pip[adjacent] != -1 and (toi, pip[adjacent]) not in self._mol.bondset:
                return False
        return True

    def _cycle_permuter(self, cycle, pip):

        def recursive_permute(pip, current_atom, remainder):
            if not remainder:
                if self._is_legal(pip, current_atom, first_atom):
                    assert (pip[current_atom] == -1)
                    pip[current_atom] = first_atom
                    yield pip
                    pip[current_atom] = -1
            else:
                for next_atom in remainder:
                    if self._is_legal(pip, current_atom, next_atom):
                        assert (pip[current_atom] == -1)
                        pip[current_atom] = next_atom
                        next_remainder = list(remainder)
                        next_remainder.remove(next_atom)
                        yield from recursive_permute(pip, next_atom, next_remainder)
                        pip[current_atom] = -1

        first_atom = cycle[0]  # First atom of the necklace
        yield from recursive_permute(pip, first_atom, cycle[1:])

    def _group_permuter(self, group, pip):
        def generate(cycles, pip):
            if not cycles:
                yield pip
            else:
                cycle = cycles[0]
                for perm in self._cycle_permuter(cycle, pip):
                    yield from generate(cycles[1:], perm)

        for cycles in self._get_cycle_structs(group):
            yield from generate(cycles, pip)

    def permute(self):
        # permutes molecule by groups
        def recursive_permute(groups, pip):
            if not groups:
                yield pip
            else:
                group = groups[0]
                groups_left = groups[1:]
                for perm in self._group_permuter(group, pip):
                    yield from recursive_permute(groups_left, perm)

        start_perm = self._empty_perm
        Groups = self._mol.equivalence_classes
        yield from recursive_permute(Groups, start_perm)
        self._num += 1


class MoleculeLegalPermuter:
    '''
    This class builds a permutation atom by atom, checking with each atom whether its new position creates an illegal permutation (one where an atom is not connected ot its neighbors)
    To that end, the class uses a list called pip (perm in progress)
    The pip is created stage by stage-- each equivalency group is built atom-by-atom (into legal cycles)
    '''

    def __init__(self, mol, opOrder, is_SN):
        self._perm_count = 0
        self._mol = mol
        self._cycle_lengths = (1, opOrder)
        if is_SN:
            self._cycle_lengths = (1, 2, opOrder)
        self._max_length = max(self._cycle_lengths)
        self._empty_perm = [-1] * (len(mol.atoms))

    def _is_legal(self, pip, to, fro):
        # fromi,j->toi,p(j)
        for adjacent in self._mol.atoms[fro].adjacent:
            if pip[adjacent] != -1 and (to, pip[adjacent]) not in self._mol.bondset:
                return False
        return True

    def _group_permuter(self, group, pip):
        '''
        permutes a single equivalency group within a molecule.
        recursively builds legal cycles within the group, atom-by-atom and cycle-by-cycle
        :return: perms-in-progress, with legal permutations of the group provided
        '''

        def do_switch(to, fro):
            if self._is_legal(pip, to, fro):
                assert pip[to] == -1
                pip[to] = fro
                return True
            return False

        def undo_switch(to, fro):
            assert pip[to] == fro
            pip[to] = -1

        def recursive_permute(pip, curr_atom, cycle_head, cycle_length, remainder):
            if cycle_length in self._cycle_lengths:  # if we have an allowed cycle length
                if do_switch(curr_atom, cycle_head):  # complete the cycle (close ends of necklace)
                    if not remainder:  # perm has been completed
                        yield pip
                    else:  # cycle has been completed, start a new cycle with remaining atoms
                        yield from recursive_permute(pip, remainder[0], remainder[0], 1, remainder[1:])
                    # roll back changes
                    undo_switch(curr_atom, cycle_head)

            # yield from permutations using a longer cycle (if applicable)
            if cycle_length < self._max_length:
                for next_atom in remainder:
                    if do_switch(curr_atom, next_atom):
                        next_remainder = list(remainder)
                        next_remainder.remove(next_atom)
                        yield from recursive_permute(pip, next_atom, cycle_head, cycle_length + 1, next_remainder)
                        # roll back changes
                        undo_switch(curr_atom, next_atom)

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

        Groups = self._mol.equivalence_classes
        yield from recursive_permute(Groups, list(self._empty_perm))


class MoleculeLegalInvertedPermuter:
    '''
    This class builds a permutation atom by atom, checking with each atom whether its new position creates an illegal permutation
    (one where an atom is not connected ot its neighbors)
    To that end, the class uses a class called pip (perm in progress)
    The pip is created stage by stage-- each equivalency group is built atom-by-atom (into legal cycles)
    '''

    class pip:
        '''
        class contains:
        permutation in progress p
        inverse permutation in progress q
        and the molecule being permuted
        '''
        def __init__(self, mol):
            self.p = [-1] * len(mol.atoms)
            self.q = [-1] * len(mol.atoms)
            self.mol = mol

        def switch(self, to, fro):
            if self._is_legal(to, fro):
                assert self.p[to] == -1 and self.q[fro] == -1
                self.p[to] = fro
                self.q[fro] = to
                return True
            return False

        def unswitch(self, to, fro):
            assert self.p[to] == fro and self.q[fro] == to
            self.p[to] = -1
            self.q[fro] = -1

        def _is_legal(self, to, fro):
            for adjacent in self.mol.atoms[fro].adjacent:
                if self.p[adjacent]!=-1 and (to, self.p[adjacent]) not in self.mol.bondset:
                    return False
            for adjacent in self.mol.atoms[to].adjacent:
                if self.q[adjacent]!=-1 and (fro, self.q[adjacent]) not in self.mol.bondset:
                    return False
            return True

    def __init__(self, mol, op_order, is_SN):
        self._perm_count = 0
        self._groups = mol.equivalence_classes
        self._pip= self.pip(mol)
        self._cycle_lengths = (1, op_order)
        if is_SN:
            self._cycle_lengths = (1, 2, op_order)
        self._max_length = op_order

    def _group_permuter(self, group, pip):
        def recursive_permute(pip, curr_atom, cycle_head, cycle_length, remainder):
            if cycle_length in self._cycle_lengths:  # if we have an allowed cycle length:
                if pip.switch(curr_atom, cycle_head):  # complete the cycle (close ends of necklace)
                    if not remainder:  # perm has been completed
                        yield pip
                    else:  # cycle has been completed, start a new cycle with remaining atoms
                        yield from recursive_permute(pip, remainder[0], remainder[0], 1, remainder[1:])
                    pip.unswitch(curr_atom, cycle_head)

            # yield from permutations using a longer cycle (if applicable)
            if cycle_length < self._max_length:
                for next_atom in remainder:
                    if pip.switch(curr_atom, next_atom):
                        next_remainder = list(remainder)
                        next_remainder.remove(next_atom)
                        yield from recursive_permute(pip, next_atom, cycle_head, cycle_length + 1, next_remainder)
                        pip.unswitch(curr_atom, next_atom)

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
            yield pip.p


class SinglePermPermuter:
    """ A permuter that returns just one permutation, used for when the permutation is specified by the user """

    def __init__(self, perm):
        self._perm = perm

    def permute(self):
        yield self._perm


if __name__ == '__main__':
    print("Testing permuters")
    # TODO: add testing code here?
