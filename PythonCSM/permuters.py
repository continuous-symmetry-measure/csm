import itertools
import colorama
from permutations import group_permuter
from permutations.utils import cycle_decomposition

__author__ = 'zmbq'

from CPP_wrapper import csm

colorama.init()

# Utility functions

# Perm generation functions
#
# perm_size: Size of permutation
# cycle_sizes: a list of legal cycle sizes
# cycle_structs: A list of cycles in a permutation: [[0], [1, 3, 4], [2], ...]

def compare(perm_size, group_size, add_groups_of_two):
    csm_perms = list(csm.GetPermuterPermutations(perm_size, group_size, add_groups_of_two))
    our_perms = list(group_permuter(perm_size, group_size, add_groups_of_two))

    allowed_cycles = {1, group_size}
    if add_groups_of_two:
        allowed_cycles.add(2)

    print("C++ permutations: %d, Python permutations: %d, unique Python: %d" % (len(csm_perms), len(our_perms), len(set(our_perms))))

    if set(csm_perms)==set(our_perms):
        return

    print("There are problems!!!")
    print("Our permutations:")
    for i, perm in enumerate(our_perms):
        cycles = cycle_decomposition(perm)
        bad_cycle = bad_perm = duplicate = False
        for cycle_size in cycles:
            if not len(cycle_size) in allowed_cycles:
                bad_cycle = True
        if not perm in csm_perms:
            bad_perm = True
        if i > 0 and perm in our_perms[:i-1]:
            duplicate = True

        print("%s %4d: %s\torder=%d\tcycles=%s" % ("Py ", i, perm, perm_order(perm), cycle_decomposition(perm)), end='')
        if bad_perm or bad_cycle:
            print(colorama.Fore.LIGHTRED_EX, end='')
            if bad_cycle:
                print('   BAD CYCLE', end='')
            elif bad_perm:
                print('   NOT IN C++', end='')
            elif duplicate:
                print('   DUPLICATE', end='')
        print(colorama.Fore.RESET)


    #print()
    #print("C++ permutations")
    for i, perm in enumerate(csm_perms):
        bad_perm = not perm in our_perms
        print("%s %4d: %s\torder=%d\tcycles=%s" % ("C++", i, perm, perm_order(perm), cycle_decomposition(perm)), end='')
        if bad_perm:
            print(colorama.Fore.LIGHTRED_EX + "   NOT IN PYTHON" + colorama.Fore.RESET, end='')
        print()


#print(list(_all_circles((2,3))))

#structs = _get_cycle_structs(5, [1,3])
#for s in structs:
#    print(s)

compare(10, 6, True)
# print(list(_all_circles((0,1,2,3))))
# print (list(_all_perms_from_cycle_struct(4, [[0,1], [2,3]])))
