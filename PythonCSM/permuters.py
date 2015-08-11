__author__ = 'zmbq'

from CPP_wrapper import csm

def apply_perm(array, perm):
    result = [0] * len(array)
    for i, p in enumerate(perm):
        # i ==> perm[i]
        result[perm[i]] = array[i]
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

def display_perms(perms):
    for i, perm in enumerate(perms):
        print("%3d: %s\torder=%d\tcycles=%s" % (i, perm, perm_order(perm), cycle_decomposition(perm)))


perms = csm.GetPermutations(8, 3, False)
display_perms(perms)

#csm.TestPermuter(5, 3, True)
#csm.TestPermuter(5, 4, True)
#csm.TestPermuter(5, 5, True)

