__author__ = 'zmbq'

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
