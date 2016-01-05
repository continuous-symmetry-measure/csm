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
    def __init__(self, _mol, opOrder, is_SN):
        self.count=0
        self.mol=_mol
        self._circle_cache = {}  # perm_size->all circles of size
        self._CACHE_LIMIT = 10
        self.cycle_size=opOrder
        self.add_cycles_of_two=is_SN
        self.constrain_set={}
    def molecule_permuter(self, elements):
        """
        Generates all permutations of a molecule
        :param elements: List of atoms being permuted- only when there are linked chains, is this is not 12345...
        :return: Generator for all the permutations
        """
        def generate(ordered_elements, perm, groups_left):
            if not groups_left: #finished permuting
                yield tuple(ordered_elements)
                return
            group = groups_left[0]
            groups_left = groups_left[1:]
            if len(group) == 1: #no options to permute for this group, continue to next
                yield from generate(ordered_elements, groups_left)
            else:
                start_elements_order=ordered_elements #this is the order the permutation was in before we called anything
                for group_perm in self.group_permuter(self, group, perm):
                    #apply perm to ordered elements
                    for i in range(len(group_perm)):
                        ordered_elements[group[i]]=start_elements_order[group_perm[i]]





        perm=np.ones(len(elements)) * -1 #permutation is empty
        yield from generate(elements, perm, self.mol.equivalence_classes)

    def group_permuter(self, group, perm):
        """
        Generates all permutations of an equivalency group
        :param group: equivalency group being permuted
        :param perm: the permutation as it has been built so far (starts as array of -1)
        :return: Generator for all the permutations
        """
        def generate():
            print("hi")
        print("hi")

    def cycle_permuter(self, cycle, perm):
        """
        Generates all permutations of a cycle (within an equivalency group)

        :param cycle: cycle being permuted
        :param perm: the permutation as it has been built so far (starts as array of -1)
        :return: Generator for all the permutations
        """
        def generate():
            print("hi")
        print("hi")

    def cycle_generator(self, group):
        """
        Generates all cycles in equivalency group of allowed size
        :param group: equivalency group cycles are being generated for
        :return: Generator for all the permutations
        """
        def generate():
            print("hi")











def ptester():
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
        elements=[i for i in range (len(mol.atoms))]
        P= permutations(mol, 2, False)
        for perm in P.molecule_permuter(elements):
            print(perm)
            NUMPERM=NUMPERM+1
        return NUMPERM


    filename=r"C:\Users\devora.witty\Sources\csm\test_cases\molDevBuilt\diamongWithArms.csm"
    num= run_test(filename)
    print(num)


    filename= r"C:\Users\devora.witty\Sources\csm\test_cases\molDevBuilt\triangleWithFroglegs.csm"
    num= run_test(filename)
    print(num)

    filename= r"C:\Users\devora.witty\Sources\csm\test_cases\molDevBuilt\triangleWithWebbedFroglegs.csm"
    num= run_test(filename)
    print(num)

    filename= r"C:\Users\devora.witty\Sources\csm\test_cases\molDevBuilt\diamondwhoseArmsHaveFingers.csm"
    num= run_test(filename)
    print(num)

    filename= r"C:\Users\devora.witty\Sources\csm\test_cases\molDevBuilt\somekindatree.csm"
    num= run_test(filename)
    print(num)

ptester()