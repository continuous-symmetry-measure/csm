import itertools

__author__ = 'zmbq'

from CPP_wrapper import csm

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

def find_cycles(perm_size, cycle_sizes):
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

def perm_from_cycles(perm_size, cycles):
    perm = list(range(perm_size))
    for cycle in cycles:
        if len(cycle)>1:
            for i in range(len(cycle)):
                perm[cycle[i]] = cycle[(i+1)%len(cycle)]
    return perm

def all_perms_in_cycle(perm, cycle_sizes):
    trivial = list(range(len(perm)))
    cycle_sizes = [len(cs) for cs in cycle_sizes if len(cs) > 1]
    i = 1
    cur = perm
    while cur != trivial:
        skip = False
        for cs in cycle_sizes:
            if (i % cs)==0:
                skip = True
                break
        if not skip:
            yield cur
        cur = apply_perm(cur, perm)
        i += 1

def get_permutations(perm_size, group_size, add_groups_of_two):
    cycle_lengths = {1, group_size}
    if add_groups_of_two:
        cycle_lengths.add(2)
    cycle_lengths = sorted(cycle_lengths)

    yield list(range(perm_size)) # The identity permutation first
    for cycle_set in find_cycles(perm_size, cycle_lengths):
        perm = perm_from_cycles(perm_size, cycle_set)
        yield from all_perms_in_cycle(perm, cycle_set)

cycles = find_cycles(6, [1, 3])
for cycle_set in cycles:
    print([cycle for cycle in cycle_set if len(cycle)>1])

perms = list(get_permutations(6, 3, False))
display_perms("Py ", perms)

print()
cpp_perms = list(csm.GetPermutations(6, 3, False))
# perms = one_cycle(5, 3)
display_perms("C++", cpp_perms)

print()
for cp in cpp_perms:
    if cp not in perms:
        print("Missing: %s" % cp)


#csm.TestPermuter(5, 3, True)
#csm.TestPermuter(5, 4, True)
#csm.TestPermuter(5, 5, True)

