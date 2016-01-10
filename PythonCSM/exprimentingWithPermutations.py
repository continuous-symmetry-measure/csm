import input_output.readers as ior
import calculations.molecule as mpy
import numpy as np
import itertools
import warnings


def first_experiment():
    #(a,b)->(p(a), p(b)
    def is_legal_permutation(perm, mol):
        m= mpy.Molecule(mol.atoms, mol.equivalence_classes)
        for atom_id in perm:
            for adj in m.atoms[atom_id].adjacent:
                if (perm[atom_id], perm[adj]) not in m.bondset:
                    return False
        return True


    def run_test(filename, count, perm):
        args_dict = {"useformat": False, "babelBond": False,
                    "inFileName": "../../test_cases/test3/c_in_1282148276_benzene.mol",
                    "ignoreSym": False, "useMass": True}

        file= open(filename, "r")

        atoms = ior.read_csm_file(file, args_dict)
        mol= mpy.Molecule(atoms)
        mol.find_equivalence_classes()

        result = is_legal_permutation(perm, mol)

        print ("done with test", count, " -- result is:", result)


    def basic_perm_concept(mol):

        #create permutation of -1
        perm = np.ones(len(mol.atoms))
        perm *= -1

        for i in range(len(mol.atoms)):
            hi= "hi"


    def permute_recursive(perm, mol, atom_index):
        curr_atom=mol.atoms[atom_index]
        for e in range(len(curr_atom.equivalency)):
            perm[atom_index]=e
            for link in range(len(curr_atom.adjacent)):
                perm[link]



    def tester():

        filename=r"C:\Users\devora.witty\Sources\csm\test_cases\molDevBuilt\diamongWithArms.csm"

        perm=[0,1,2,3,4,5]
        run_test(filename, 1, perm)
        perm=[1,2,3,4,5,0]
        run_test(filename, 2, perm)

        filename= r"C:\Users\devora.witty\Sources\csm\test_cases\molDevBuilt\triangleWithFroglegs.csm"
        perm=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
        run_test(filename, 3, perm)
        perm=[8,9,10,11, 0,1,2,3, 4,5,6,7,14,12,13,17,15,16]
        run_test(filename, 4, perm)

        print("done")


class old_permutations:
    def __init__(self, _mol, opOrder, is_SN):
        self.count=0
        self.mol=_mol
        self._circle_cache = {}  # perm_size->all circles of size
        self._CACHE_LIMIT = 10
        self.cycle_size=opOrder
        self.add_cycles_of_two=is_SN
        self.constrain_set={}

    def is_legal_atom_permutation(self, perm, atom_id):
        for adj in self.mol.atoms[atom_id].adjacent:
            if (perm[atom_id], perm[adj]) not in self.mol.bondset and perm[adj]!=-1:
                return False
        return True

    def is_legal_cycle_permutation(self, perm, cycle):
        for atom_id in cycle:
            for adj in self.mol.atoms[atom_id].adjacent:
                if (perm[atom_id], perm[adj]) not in self.mol.bondset:
                    return False
        return True

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

                if cycle_size==1:
                    yield from generate(cycles, rest)
                    continue

                start = (left[0],)
                for rest in itertools.combinations(rest, cycle_size-1):
                    cycle = start + rest
                    new_left = [item for item in left if item not in cycle]
                    yield from generate(cycles+[cycle], new_left)

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
        trivial = tuple(range(1, size))   # (1,2,...size-1)
        for necklace in itertools.permutations(trivial):
            # The actual necklace is [0]+necklace
            cycle_perm = [0] * size
            cur = 0
            for element in necklace:
                cycle_perm[cur] = element
                cur = element
            cycle_perm[necklace[-1]] = 0  # Add the [0] that is missing from the necklace
            yield cycle_perm

    def all_circle_permutations(self, size):
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
                yield from generate(perm, cycle_index-1)
            else:
                # Example:
                # Lets say the permutation is (0, 1 ,2 ,3), and the cycle is (0, 1, 3)
                circles = self.all_circle_permutations(len(cycle))
                for circle_perm in circles:
                    # _all_circle_permtuations yields (1, 2, 0) and (2, 0 ,1)
                    # The permutations we need to return are (1, 3, 2, 0) and (3, 0, 2, 1) - these have
                    # one stationary point - 2, and a cycle of length 3.
                    # To do this, we need to convert the cycle from (1,2,0) and (2,0,1) to (1,3,0) and (3,0,1)
                    for i in range(len(cycle)):
                        converted_circle = cycle[circle_perm[i]]
                        perm[cycle[i]] = converted_circle  # Perm is now (1, 3, 2, 0) or (3, 0, 2, 1) <--- The converted_circle applied to (0, 1, 3)
                    yield from generate(perm, cycle_index-1)    # Apply the rest of the circles

        yield from generate(trivial[:], len(cycle_struct)-1)

    def group_permuter(self, group_size):
        """
        Generates all permutations of size and cycle sizes.
        :param group_size: Size of group
        :param cycle_size: Cycle size
        :param add_cycles_of_two: When true, cycles of size of two are legal
        :return: Generator for all the permutations
        """
        cycle_lengths = {1, self.cycle_size}
        if self.add_cycles_of_two:
            cycle_lengths.add(2)
        cycle_lengths = sorted(cycle_lengths)

        for cycle_struct in self._get_cycle_structs(group_size, cycle_lengths):
            # Loop over all cycle structures, which are of the allowed lengths covering the entire permutation
            yield from self._all_perms_from_cycle_struct(group_size, cycle_struct) # Return all permutations for each cycle structure

    def molecule_permuter(self, elements):
        """
        Generates all permutations of a molecule
        :param elements: List of atoms being permuted- when there are linked chains, this is not 12345...
        :return: Generator for all the permutations
        """
        groups=self.mol.equivalence_classes
        print("entered molecule permuter")

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
                for group_perm in self.group_permuter(len(group)):
                    # group_permuter yields (1, 2, 0) and (2, 0 ,1)
                    # The permutations we need to return are (1, 3, 2, 0) and (3, 0, 2, 1) - these have
                    # one stationary point - 2, and a cycle of length 3.
                    # To do this, we need to convert the group_perm from (1,2,0) and (2,0,1) to (1,3,0) and (3,0,1)
                    converted_group = [group[i] for i in group_perm]  # Converted circle is now (1,3,0) or (3,0,1)
                    for i in range(len(converted_group)):
                        ordered_elements[group[i]] = start_elements_order[converted_group[i]]
                        # Perm is now (1, 3, 2, 0) or (3, 0, 2, 1) <--- The converted_circle applied to (0, 1, 3)

                    yield from generate(ordered_elements, groups_left)    # Apply the rest of the circles
                    ordered_elements = start_elements_order[:]

        elements_order = elements  # The starting elements order
        yield from generate(elements_order, groups)



class permutation:
    def __init__(self, mol, opOrder=2, is_SN=False):
        self.count=0
        self._mol=mol
        self._len=len(mol.atoms)
        self._identityPerm=np.arange(self._len)
        self._is_SN=is_SN
        self._opOrder= opOrder


    def is_legal_perm(self, PiP, position, index):
        for adj in self._mol.atoms[index].adjacent:
            if (self._identityPerm[position], self._identityPerm[adj]) not in self._mol.bondset and PiP[adj]!=-1:
                return False
        return True


    def cycle_permuter(self, group, PIP):

        Indexes=list(group)
        Positions=list(Indexes)

        def generate(indexes, positions, pip):
            if positions==[]:
                yield pip, indexes
            else:
                position=positions.pop()
                for index in indexes:
                    if self.is_legal_perm(pip, position,index):
                        carryPip=np.array(pip)
                        carryPip[position]=self._identityPerm[index]
                        indexes.remove(index)
                        yield from generate(indexes, positions, carryPip)

        First_pos= Positions.pop()
        for perm, remainder in generate(Indexes, Positions, PIP):
            index=remainder.pop()
            if self.is_legal_perm(perm, First_pos, index):
                perm[First_pos]=index
                yield perm

    def cycle_creator(self, group):
        cycle_sizes={self._opOrder}
        if self._is_SN:
            cycle_sizes.add(2)

        def generate(cycles, left):
            len_left = len(left)
            if len_left == 0:
                yield cycles

            for cycle_size in cycle_sizes:
                if cycle_size > len_left:
                    continue
                rest = left[1:]

                start = (left[0],)
                for rest in itertools.combinations(rest, cycle_size-1):
                    cycle = start + rest
                    new_left = [item for item in left if item not in cycle]
                    yield from generate(cycles+[cycle], new_left)

        yield from generate([], list(group))



    def molecule_permuter(self):
        Blank_perm=np.ones(len(self._mol.atoms)) * -1
        groups=self._mol.equivalence_classes

        def generate(PIP, groups_left):
            if not groups_left:
                yield PIP

            group=groups_left[0]
            groups_left=groups_left[1:]
            for cycles in self.cycle_creator(group):
                for cycle in cycles:
                    for perm in self.cycle_permuter(cycle, PIP):
                        PIP=perm
                        groups_left=groups[1:]
                        yield from (PIP, groups_left)

        yield from generate(Blank_perm, groups)

    def simple_generator(self):
        for Perm in self.molecule_permuter():
            print (Perm)

def test():
        args_dict = {"useformat": False, "babelBond": False,
                    "inFileName": "../../test_cases/test3/c_in_1282148276_benzene.mol",
                    "ignoreSym": False, "useMass": True}
        def run_test(filename):
            print("running test")
            NUMPERM=0
            file= open(filename, "r")
            atoms = ior.read_csm_file(file, args_dict)
            mol= mpy.Molecule(atoms)
            mol.find_equivalence_classes()
            p=permutation(mol)
            p.simple_generator()

        filename=r"C:\Users\devora.witty\Sources\csm\test_cases\molDevBuilt\diamongWithArms.csm"
        run_test(filename)

class molecule_permuter:
    def __init__(self, mol, opOrder, is_SN):
        self.num=0
        self.mol=mol
        self.opOrder=opOrder
        self.add_cycles_of_two=is_SN
        self.empty_perm=[-1]*(len(mol.atoms))
        self.hi="ji"

    def _get_cycle_structs(self, group):
        cycle_sizes={1, self.opOrder}
        if self.add_cycles_of_two:
            cycle_sizes.add(2)
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

                #if cycle_size==1:
                #    yield group
                #    continue

                start = (left[0],)
                for rest in itertools.combinations(rest, cycle_size-1):
                    cycle = start + rest
                    new_left = [item for item in left if item not in cycle]
                    yield from generate(cycles+[cycle], new_left)

        return generate([], list(group))

    def is_legal(self, Pip, toi, fromi):
        #fromi,j->toi,p(j)
        for adjacent in self.mol.atoms[fromi].adjacent:
            if  Pip[adjacent]!=-1 and (toi, Pip[adjacent]) not in self.mol.bondset:
                return False
        return True

    def cycle_permuter(self, cycle, pip):

        def recursive_permute(pip, current_atom, remainder):
            if not remainder:
                if self.is_legal(pip, current_atom, first_atom):
                    # Make sure that pip[current_atom] is -1
                    pip[current_atom] = first_atom
                    yield pip
                    # Set pip[current_atom] back to -1
            else:
                for next_atom in remainder:
                    if self.is_legal(pip, current_atom, next_atom):
                        # Make sure that pip[current_atom] is -1
                        pip[current_atom] = next_atom
                        next_remainder = list(remainder)
                        next_remainder.remove(next_atom)
                        yield from recursive_permute(pip, next_atom, next_remainder)
                        # Set pip[current_atom] back to -1

        first_atom = cycle[0]  # First atom of the necklace
        yield from recursive_permute(pip, first_atom, cycle[1:])



    def group_permuter(self, group, Pip):
        def generate(cycles, Pip):
            if not cycles:
                yield Pip
            else:
                cycle=cycles[0]
                for perm in self.cycle_permuter(cycle, Pip):
                    yield from generate(cycles[1:], perm)

        for cycles in self._get_cycle_structs(group):
            yield from generate(cycles, Pip)



    def molecule_permuter(self):
        #permutes molecule by groups
        def recursive_permute(groups, Pip):
            if not groups:
                yield Pip
            else:
                group=groups[0]
                groups_left=groups[1:]
                for perm in self.group_permuter(group, Pip):
                    yield from recursive_permute(groups_left, perm)

        Pip=self.empty_perm
        Groups= self.mol.equivalence_classes
        yield from recursive_permute(Groups, Pip)
        self.num+=1







def run_tests():
        args_dict = {"useformat": False, "babelBond": False,
                    "inFileName": "../../test_cases/test3/c_in_1282148276_benzene.mol",
                    "ignoreSym": False, "useMass": True}
        def run_test(filename):
            print("running test", filename)
            NUMPERM=0
            file= open(filename, "r")
            atoms = ior.read_csm_file(file, args_dict)
            mol= mpy.Molecule(atoms)
            mol.find_equivalence_classes()
            p= molecule_permuter(mol, 3, True)
            elements=[i for i in range (len(mol.atoms))]
            Pip=np.ones(len(elements))*-1
            for perm in p.molecule_permuter():
                print(perm)
                NUMPERM +=1
            return NUMPERM

        filename=r"C:\Users\devora.witty\Sources\csm\test_cases\molDevBuilt\diamongWithArms.csm"
        num= run_test(filename)
        print(num)

        filename= r"C:\Users\devora.witty\Sources\csm\test_cases\molDevBuilt\somekindatree.csm"
        num= run_test(filename)
        print(num)
