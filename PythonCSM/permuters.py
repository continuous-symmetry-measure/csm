import itertools
import sys
from timeit import Timer
from calculations.csm_calculations_data import CSMCalculationsData
from calculations.preprocess_molecule import preprocess_molecule
import colorama
from permutations import group_permuter, molecule_permuter
from permutations.permuters import _len_group_permuter, _get_cycle_structs, _all_circle_permutations
from permutations.utils import cycle_decomposition, perm_order
from arguments import process_arguments, create_parser
import numpy as np

__author__ = 'zmbq'

from CPP_wrapper import csm, permutations, experiments

colorama.init()


# Utility functions


# Perm generation functions
#
# perm_size: Size of permutation
# cycle_sizes: a list of legal cycle sizes
# cycle_structs: A list of cycles in a permutation: [[0], [1, 3, 4], [2], ...]
def compare(perm_size, group_size, add_groups_of_two):
    csm_perms = []
    for csm_perm in csm.GetPermuterPermutations(perm_size, group_size, add_groups_of_two):
        csm_perms.append(tuple(csm_perm))

    our_perms = []
    for our_perm in permutations.group_permuter(perm_size, group_size, add_groups_of_two):
        our_perms.append(tuple(our_perm))

    allowed_cycles = {1, group_size}
    if add_groups_of_two:
        allowed_cycles.add(2)

    print("C++ permutations: %d, Python permutations: %d, unique Python: %d" % (
    len(csm_perms), len(our_perms), len(set(our_perms))))

    csm_perms = [tuple(perm) for perm in csm_perms]
    if set(csm_perms) == set(our_perms):
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
        if i > 0 and perm in our_perms[:i - 1]:
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


    # print()
    # print("C++ permutations")
    for i, perm in enumerate(csm_perms):
        bad_perm = perm not in our_perms
        print("%s %4d: %s\torder=%d\tcycles=%s" % ("C++", i, perm, perm_order(perm), cycle_decomposition(perm)), end='')
        if bad_perm:
            print(colorama.Fore.LIGHTRED_EX + "   NOT IN PYTHON" + colorama.Fore.RESET, end='')
        print()


def compare_molecule(args):
    parser = create_parser()
    result = parser.parse_args(args)
    csm_args = process_arguments(result)
    preprocess_molecule(csm_args)
    csm.SetCSMOptions(csm_args)

    csm_perms = list(csm.GetMoleculePermutations())
    our_perms = list(csm.molecule_permuter(len(csm_args['molecule'].atoms), csm_args['molecule'].equivalence_classes,
                                       csm_args['opOrder'], csm_args['type'] == 'SN'))

    print("C++ permutations: %d, Python permutations: %d, unique Python: %d" % (
    len(csm_perms), len(our_perms), len(set(our_perms))))

    if set(csm_perms) == set(our_perms):
        print("OK")
    else:
        print("There are problems!!!")
    """
        print("Converting ...")

        convert_perm = []
        for i in range(len(our_perms[0])):
            convert_perm.append(csm_perms[len(csm_perms)-1][our_perms[0][i]])

        our_perms = [list(perm) for perm in our_perms]
        for j in range(len(our_perms)):
            for i in range(len(our_perms[j])):
                our_perms[j][i] = convert_perm[our_perms[j][i]]

        our_perms = [tuple(perm) for perm in our_perms]
        if set(csm_perms) == set(our_perms):
            print("Now it's OK")
        else:
            print("There are problems still!!!")
    """
    print("Groups:")
    for group in csm_args['molecule'].equivalence_classes:
        print("\t%s" % group)

    print("Our permutations:")
    for i, perm in enumerate(our_perms):
        bad_perm = duplicate = False
        if not perm in csm_perms:
            bad_perm = True
        if perm in our_perms[:i-1]:
            duplicate = True

        print("%s %4d: %s" % ("Py ", i, perm), end='')
        if bad_perm:
            print(colorama.Fore.LIGHTRED_EX, end='')
            if bad_perm:
                print('   NOT IN C++', end='')
            elif duplicate:
                print('   DUPLICATE', end='')
        print(colorama.Fore.RESET)


    # print()
    # print("C++ permutations")
    for i, perm in enumerate(csm_perms):
        bad_perm = perm not in our_perms
        print("%s %4d: %s" % ("C++", i, perm), end='')
        if bad_perm:
            print(colorama.Fore.LIGHTRED_EX + "   NOT IN PYTHON" + colorama.Fore.RESET, end='')
        print()

    data = CSMCalculationsData(csm_args)
    result = csm.CsmOperation(data)

def big_test():
    group_size, cycle_size, add_groups_of_two = 12, 8, True
    cycle_sizes = {1, cycle_size}
    if add_groups_of_two:
        cycle_sizes.add(2)

    def count_cpp():
        count = 0
        for perm in csm.GetPermuterPermutations(group_size, cycle_size, add_groups_of_two):
            count += 1
        print('C++ count: %d' % count)

    def count_python():
        count = 0
        #for struct in _get_cycle_structs(group_size, cycle_sizes):
        for perm in group_permuter(group_size, cycle_size, add_groups_of_two):
            count += 1
        print('Python count: %d' % count)


    def count_cython():
        count = 0
        #for struct in permutations._get_cycle_structs(group_size, cycle_sizes):
        for perm in permutations.group_permuter(group_size, cycle_size, add_groups_of_two):
            count += 1
        print('Cython count: %d' % count)

    timer_cython = Timer(count_cython)
    time_cython = timer_cython.timeit(number=1)
    print("Cython: %s" % time_cython)

    #timer_python = Timer(count_python)
    #time_python = timer_python.timeit(number=1)
    #print("Python: %s" % time_python)

#    timer_cpp = Timer(count_cpp)
#    time_cpp = timer_cpp.timeit(number=1)
#    print("C++: %s" % time_cpp)

big_test()

# compare(10, 5, True)

def some_experiment():
    def experiment(func, size, count):
        total = 0
        for entry in func(size, count):
            # print(entry)
            total += 1

    funcs = [#experiments.measure_straight,
             #experiments.measure_typed,
             experiments.measure_array,
             experiments.measure_memory_view,
             experiments.measure_pointer,
             experiments.measure_numpy]

    for func in funcs:
        timer = Timer(lambda: experiment(func, 5, 4000000))
        print(func.__name__, '...')
        print(timer.timeit(number=3) / 3)
#some_experiment()

# print(list(_all_circles((2,3))))

# structs = _get_cycle_structs(5, [1,3])
# for s in structs:
#    print(s)

def count_structs(perm_size, cycle_sizes):
    count = 0
    for struct in _get_cycle_structs(perm_size, cycle_sizes):
        count += 1
    return count

def ratio(group_size, cycle_size, add_groups_of_two):
    cycle_sizes = {1, cycle_size}
    if add_groups_of_two:
        cycle_sizes.add(2)
    structs = count_structs(group_size, cycle_sizes)
    total = _len_group_permuter(group_size, cycle_size, add_groups_of_two)
    return total / structs

#print(ratio(12, 2, True))
compare(11, 6, True)
# print(list(_all_circles((0,1,2,3))))
# print (list(_all_perms_from_cycle_struct(4, [[0,1], [2,3]])))

# compare_molecule(sys.argv[1:])
