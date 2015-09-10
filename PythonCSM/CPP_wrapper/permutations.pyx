__author__ = 'zmbq'

import itertools
from cpython cimport array
cimport cython

cdef array.array _int_array_template = array.array('i', [])

def _get_cycle_structs(perm_size, cycle_sizes):
    """
    Generates a list of cycles in a permutation. The cycles cover the entire permutation, and are only of sizes in cycle_sizes
    :param perm_size: Permutation size
    :param cycle_sizes: Allowed cycle sizes
    :return: A generator of a list of cycles covering the entire permtutation

    """
    def generate(cycles, left):
        cdef array.array cycle
        cdef int[:] cycle_buf
        cdef int i

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
                cycle = array.clone(_int_array_template, cycle_size, False)
                # cycle_buf = cycle
                cycle[0] = left[0]
                for i in range(1, cycle_size):
                    cycle[i] = rest[i-1]


                new_left = [item for item in left if item not in cycle]
                yield from generate(cycles+[cycle], new_left)

    if 1 not in cycle_sizes:
        raise ValueError("1 must be a valid cycle size")
    yield from generate([], list(range(perm_size)))


def _all_circle_permutations(size):
    """ Returns the permutation of the cycle """
    # To compute a full cycle of length n, we take a permutation p of size n-1 and
    # create the cycle like so: 0 goes to p[0], p[0] to p[1], p[1] to p[2] and so on, until p[n-1] goes back to 0
    # For this to work, p needs to be a permutation of 1..n-1.
    #
    # For example, size = 3 we have p=(1,2) for the circle 0->1->2->0 and p=(2,1) for the circle 0->2->1->0
    cdef array.array cycle_perm = array.clone(_int_array_template, size, False)
    cdef int[] cycle_perm_ptr = cycle_perm.data.as_ints
    cdef int element, cur

    trivial = tuple(range(1, size))   # (1,2,...size-1)
    for necklace in itertools.permutations(trivial):
        # The actual necklace is [0]+necklace
        cur = 0
        for element in necklace:
            cycle_perm_ptr[cur] = element
            cur = element
        cycle_perm_ptr[cur] = 0  # Add the [0] that is missing from the necklace
        yield cycle_perm_ptr



def _all_perms_from_cycle_struct(perm_size, cycle_struct):
    """
    Generates all permutations with a specific cycle structure
    :param perm_size: Permutation size
    :param cycle_struct: Cycle structure (list of cycles, each cycle is itself a list of its indices)
    :return: Generator that generates all the permutations
    """
    cdef array.array perm = array.clone(_int_array_template, perm_size, False)
    cdef int[:] perm_buf = perm
    cdef array.array temp = array.clone(_int_array_template, perm_size, False)
    cdef int[:] temp_buf = temp
    cdef int i

    for i in range(perm_size):
        perm_buf[i] = i

    def generate(int cycle_index):
        # Goes over all the circles of the first cycle, apply each circle and
        # recursively generates the circles of the rest of the cycles
        # perm is built gradually, each recursive call applies one cycle
        cdef int i
        # cdef int[:] cycle, circle_perm

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
            temp_buf[:] = perm_buf
            for circle_perm in _all_circle_permutations(cycle_len):
                # _all_circle_permtuations yields (1, 2, 0) and (2, 0 ,1)
                # The permutations we need to return are (1, 3, 2, 0) and (3, 0, 2, 1) - these have
                # one stationary point - 2, and a cycle of length 3.
                # To do this, we need to convert the cycle from (1,2,0) and (2,0,1) to (1,3,0) and (3,0,1)
                converted_circle = [cycle[i] for i in circle_perm]  # Converted circle is now (1,3,0) or (3,0,1)
                for i in range(len(converted_circle)):
                    perm_buf[cycle[i]] = converted_circle[i]  # Perm is now (1, 3, 2, 0) or (3, 0, 2, 1) <--- The converted_circle applied to (0, 1, 3)
                yield from generate(cycle_index-1)    # Apply the rest of the circles
                perm_buf[:] = temp_buf

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

