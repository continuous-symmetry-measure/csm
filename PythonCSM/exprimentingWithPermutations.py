import input_output.readers as ior
import calculations.molecule as mpy
import numpy as np
import itertools
import warnings


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
                    assert(pip[current_atom]==-1)
                    pip[current_atom] = first_atom
                    yield pip
                    pip[current_atom]=-1
            else:
                if pip[current_atom]==-1:
                    for next_atom in remainder:
                        if self.is_legal(pip, current_atom, next_atom):
                            assert(pip[current_atom]==-1)
                            pip[current_atom] = next_atom
                            next_remainder = list(remainder)
                            next_remainder.remove(next_atom)
                            yield from recursive_permute(pip, next_atom, next_remainder)
                            pip[current_atom]= -1

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
            p= molecule_permuter(mol, 3, False)
            elements=[i for i in range (len(mol.atoms))]
            Pip=np.ones(len(elements))*-1
            for perm in p.molecule_permuter():
                NUMPERM +=1
            return NUMPERM


        filename=r"C:\Users\devora.witty\Sources\csm\test_cases\test4\ZnCu10-p9.csm"
        num= run_test(filename)
        print(num)

        filename= r"C:\Users\devora.witty\Sources\csm\test_cases\molDevBuilt\somekindatree.csm"
        num= run_test(filename)
        print(num)

run_tests()