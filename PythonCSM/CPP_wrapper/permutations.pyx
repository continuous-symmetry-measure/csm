__author__ = 'zmbq'

import itertools
from cpython cimport array
cimport cython

cdef array.array _int_array_template = array.array('i', [])
cdef inline int *ptr(array.array int_array):
    return int_array.data.as_ints

cdef array.array int_array(int size, zeros=False):
    return array.clone(_int_array_template, size, zeros)

def _get_cycle_structs(perm_size, cycle_sizes):
    """
    Generates a list of cycles in a permutation. The cycles cover the entire permutation, and are only of sizes in cycle_sizes
    :param perm_size: Permutation size
    :param cycle_sizes: Allowed cycle sizes
    :return: A generator of a list of cycles covering the entire permtutation

    """
    def generate(cycles, left):
        cdef array.array cycle
        cdef int i
        cdef int *p_cycle

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

            for rest in itertools.combinations(rest, cycle_size-1):
                cycle = int_array(cycle_size)
                p_cycle = ptr(cycle)

                p_cycle[0] = left[0]
                for i in range(1, cycle_size):
                    p_cycle[i] = rest[i-1]

                new_left = [item for item in left if item not in cycle]
                yield from generate(cycles+[cycle], new_left)

    if 1 not in cycle_sizes:
        raise ValueError("1 must be a valid cycle size")
    yield from generate([], list(range(perm_size)))


def _calc_all_circle_permutations(size):
    """ Returns the permutation of the cycle """
    # To compute a full cycle of length n, we take a permutation p of size n-1 and
    # create the cycle like so: 0 goes to p[0], p[0] to p[1], p[1] to p[2] and so on, until p[n-1] goes back to 0
    # For this to work, p needs to be a permutation of 1..n-1.
    #
    # For example, size = 3 we have p=(1,2) for the circle 0->1->2->0 and p=(2,1) for the circle 0->2->1->0
    cdef array.array cycle_perm = int_array(size)
    cdef int *p_cycle_perm = ptr(cycle_perm)
    cdef int element, cur

    trivial = tuple(range(1, size))   # (1,2,...size-1)
    for necklace in itertools.permutations(trivial):
        # The actual necklace is [0]+necklace
        cur = 0
        for element in necklace:
            p_cycle_perm[cur] = element
            cur = element
        p_cycle_perm[cur] = 0
        yield cycle_perm

_circle_cache = {}  # perm_size->all circles of size
_CACHE_LIMIT = 10

def _all_circle_permutations(size):
    if size > _CACHE_LIMIT:
        return _calc_all_circle_permutations(size)

    if not size in _circle_cache:
        entries = list(_calc_all_circle_permutations(size))
        _circle_cache[size] = entries
    result = _circle_cache[size]
    return result


def _all_perms_from_cycle_struct(perm_size, cycle_struct):
    """
    Generates all permutations with a specific cycle structure
    :param perm_size: Permutation size
    :param cycle_struct: Cycle structure (list of cycles, each cycle is itself a list of its indices)
    :return: Generator that generates all the permutations
    """
    cdef array.array perm = int_array(perm_size)
    cdef int *p_perm = ptr(perm)

    cdef array.array orig_perm = int_array(perm_size)
    cdef int *p_orig_perm = ptr(orig_perm)

    cdef array.array converted_circle = int_array(perm_size)
    cdef int *p_converted_circle = ptr(converted_circle)

    cdef int i

    for i in range(perm_size):
        p_perm[i] = i

    def generate(int cycle_index):
        # Goes over all the circles of the first cycle, apply each circle and
        # recursively generates the circles of the rest of the cycles
        # perm is built gradually, each recursive call applies one cycle
        cdef int i
        cdef int *p_circle_perm
        cdef int *p_cycle

        if cycle_index < 0:
            yield perm
            return

        cycle = cycle_struct[cycle_index]
        cdef int cycle_len = len(cycle)
        if cycle_len == 1:
            yield from generate(cycle_index-1)
        else:
            # Example:
            # Lets say the permutation is (0, 1 ,2 ,3), and the cycle is (0, 1, 3)

            for circle_perm in _all_circle_permutations(cycle_len):
                p_circle_perm = ptr(circle_perm)
                p_cycle = ptr(cycle)
                # _all_circle_permtuations yields (1, 2, 0) and (2, 0 ,1)
                # The permutations we need to return are (1, 3, 2, 0) and (3, 0, 2, 1) - these have
                # one stationary point - 2, and a cycle of length 3.
                # To do this, we need to convert the cycle from (1,2,0) and (2,0,1) to (1,3,0) and (3,0,1)
                for i in range(cycle_len):
                    p_converted_circle[i] = p_cycle[p_circle_perm[i]]
                for i in range(cycle_len):
                    p_orig_perm[i] = p_perm[p_cycle[i]]
                    p_perm[p_cycle[i]] = p_converted_circle[i]
                yield from generate(cycle_index-1)    # Apply the rest of the circles

                for i in range(cycle_len):
                    p_perm[p_cycle[i]] = p_orig_perm[i]

    yield from generate(len(cycle_struct)-1)



def group_permuter(int group_size, int cycle_size, int add_cycles_of_two):
    """
    Generates all permutations of size and cycle sizes.
    :param group_size: Size of group
    :param cycle_size: Cycle size
    :param add_cycles_of_two: When true, cycles of size of two are legal
    :return: Generator for all the permutations
    """
    cycle_lengths = {1, cycle_size}
    if add_cycles_of_two:
        cycle_lengths.add(2)
    cycle_lengths = sorted(cycle_lengths)

    for cycle_struct in _get_cycle_structs(group_size, cycle_lengths):
        # Loop over all cycle structures, which are of the allowed lengths covering the entire permutation
        yield from _all_perms_from_cycle_struct(group_size, cycle_struct) # Return all permutations for each cycle structure


def molecule_permuter(int molecule_size, int groups, int cycle_size, int add_cycles_of_two):
    """
    Generates all permutations of a molecule
    :param molecule_size: Molecule size
    :param groups: Equivalency groups
    :param cycle_size: Allowed cycle size
    :param add_cycles_of_two: When true, cycles of size of two are legal
    :return: Generator for all the permutations
    """
    def generate(int[:]perm_buf, groups_left):
        """ Goes over all the permutations of the first group, applies each
         permutation and recursively aplies permutations of the entire groups
        :param perm:
        :param groups_left:
        :return:
        """
        cdef int[:] group_perm

        if not groups_left:
            yield perm_buf
            return

        group = groups_left[0]
        groups_left = groups_left[1:]
        if len(group) == 1:
            yield from generate(perm_buf, groups_left)
        else:
            # Example:
            # Lets say the permutation is (0, 1 ,2 ,3), and the group is (0, 1, 3)
            temp_buf[:] = perm
            for group_perm in group_permuter(len(group), cycle_size, add_cycles_of_two):
                # group_permuter yields (1, 2, 0) and (2, 0 ,1)
                # The permutations we need to return are (1, 3, 2, 0) and (3, 0, 2, 1) - these have
                # one stationary point - 2, and a cycle of length 3.
                # To do this, we need to convert the group_perm from (1,2,0) and (2,0,1) to (1,3,0) and (3,0,1)
                converted_group = [group[i] for i in group_perm]  # Converted circle is now (1,3,0) or (3,0,1)
                for i in range(len(converted_group)):
                    perm[group[i]] = converted_group[i]  # Perm is now (1, 3, 2, 0) or (3, 0, 2, 1) <--- The converted_circle applied to (0, 1, 3)
                yield from generate(perm, groups_left)    # Apply the rest of the circles
                perm[:] = temp_buf

    perm = array.array('i', list(range(molecule_size)))  # The trivial permutation
    temp = array.clone(perm, molecule_size, False)
    cdef int[:] perm_buf = perm
    cdef int[:] temp_buf = temp
    yield from generate(perm_buf, groups)

