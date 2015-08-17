import itertools

__author__ = 'zmbq'

from CPP_wrapper import csm

# Utility functions
def apply_perm(array, perm):
    result = [0] * len(array)
    for i, p in enumerate(perm):
        result[i] = array[perm[i]]
    return result

def perm_order(perm):
    start = list(range(-len(perm), 0))
    p = apply_perm(start, perm)
    order = 1
    while p!=start:
        p = apply_perm(p, perm)
        order += 1
    return order

def cycle_decomposition(perm):
    unvisited = set(range(len(perm)))  # All indices are unvisited
    def find_cycle(start):
        cycle = [start]
        next = perm[start]
        while next!=start:
            cycle.append(next)
            unvisited.remove(next)
            next = perm[next]
        return cycle

    cycles = []
    while unvisited:
        start = unvisited.pop()
        cycle = find_cycle(start)
        if len(cycle) > 1:
            cycles.append(cycle)

    return cycles

def display_perms(prefix, perms):
    for i, perm in enumerate(perms):
        print("%s %4d: %s\torder=%d\tcycles=%s" % (prefix, i, perm, perm_order(perm), cycle_decomposition(perm)))

# Perm generation functions
#
# perm_size: Size of permutation
# cycle_sizes: a list of legal cycle sizes
# cycle_structs: A list of cycles in a permutation: [[0], [1, 3, 4], [2], ...]

def get_cycle_structs(perm_size, cycle_sizes):
    def generate(cycles, left):
        len_left = len(left)
        if len_left == 0:
            yield cycles
            return

        for cycle_size in cycle_sizes:
            if cycle_size > len_left:
                break
            start = [left[0]]
            for rest in itertools.combinations(left[1:], cycle_size-1):
                cycle = start + list(rest)
                new_left = [item for item in left if item not in cycle]
                yield from generate(cycles+[cycle], new_left)

    yield from generate([], list(range(perm_size)))


def all_circular_permutations(elements):
    first = elements[0]
    rest = elements[1:]
    for rest_perm in itertools.permutations(rest):
        yield [first] + list(rest_perm)

def perm_from_cycle_struct(perm_size, cycle_structs):
    # Each cycle_struct yields one permutation - one that rotates all the cycles once.
    perm = list(range(perm_size))
    for cycle in cycle_structs:
        if len(cycle)>1:
            for i in range(len(cycle)):
                perm[cycle[i]] = cycle[(i+1)%len(cycle)]
    return perm


def all_perms_from_cycle_struct(start_point, perm, cycle_sizes):
    trivial = list(range(len(perm)))
    cycle_sizes = [len(cs) for cs in cycle_sizes if len(cs) > 1]
    i = 1

    cur = start_point
    if cur != trivial:
        yield cur
    cur = apply_perm(cur, perm)

    while cur != start_point:
        skip = False
        for cs in cycle_sizes:
            if (i % cs) == 0:
                skip = True
                break
        if not skip:
            print("%s, %s, %s, %2d: %s" % (start_point, perm, cycle_sizes, i, cur))
            yield cur
        cur = apply_perm(cur, perm)
        i += 1


def cycle_struct_start_points(perm_size, cycle_struct):
    print("cycle_struct_start_points %s, %s" % (perm_size, cycle_struct))
    trivial = list(range(perm_size))
    yield trivial
    print("\t==> %s" % trivial)

    def generate(perm, cycle_struct):
        print("cycle_struct_start_points.generate %s, %s" % (perm, cycle_struct))
        if not cycle_struct and perm!=trivial:
            print("\t==> %s" % perm)
            yield perm

        if not cycle_struct:
            return

        cycle = cycle_struct[0]
        for cycle_perm in all_circular_permutations(cycle):
            new_perm = perm[:]
            for i in range(len(cycle)):
                new_perm[cycle[i]] = cycle_perm[i]
            yield from generate(new_perm, cycle_struct[1:])

    trivial = list(range(perm_size))
    yield from generate(trivial, cycle_struct)


def get_permutations(perm_size, group_size, add_groups_of_two):
    cycle_lengths = {1, group_size}
    if add_groups_of_two:
        cycle_lengths.add(2)
    cycle_lengths = sorted(cycle_lengths)

    yield list(range(perm_size)) # The identity permutation first
    for cycle_struct in get_cycle_structs(perm_size, cycle_lengths):
        perm = perm_from_cycle_struct(perm_size, cycle_struct)
        start_points = cycle_struct_start_points(perm_size, cycle_struct)
        for start_point in start_points:
            yield from all_perms_from_cycle_struct(start_point, perm, cycle_struct)


def compare(perm_size, group_size, add_groups_of_two):
    csm_perms = list(csm.GetPermutations(perm_size, group_size, add_groups_of_two))
    our_perms = list(get_permutations(perm_size, group_size, add_groups_of_two))

    print("Our permutations:")
    display_perms("Py", our_perms)

    print()
    print("C++ permutations")
    display_perms("C++", csm_perms)

compare(4, 3, False)

#cycles = find_cycles(6, [1, 3])
#for cycle_set in cycles:
#    print([cycle for cycle in cycle_set if len(cycle)>1])



#perms = get_permutations(6, 3, False))
#display_perms("Py ", perms)

#print()
#cpp_perms = list(csm.GetPermutations(6, 3, False))
# perms = one_cycle(5, 3)
#display_perms("C++", cpp_perms)

#print()
#for cp in cpp_perms:
#    if cp not in perms:
#        print("Missing: %s" % cp)


#csm.TestPermuter(5, 3, True)
#csm.TestPermuter(5, 4, True)
#csm.TestPermuter(5, 5, True)

