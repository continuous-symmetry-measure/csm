__author__ = 'zmbq'

include "misc.pxi"
import itertools
from libcpp cimport bool

cimport cython


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


def _calc_all_circle_permutations(size, bool for_caching=False):
    """ Returns the permutation of the cycle """
    # To compute a full cycle of length n, we take a permutation p of size n-1 and
    # create the cycle like so: 0 goes to p[0], p[0] to p[1], p[1] to p[2] and so on, until p[n-1] goes back to 0
    # For this to work, p needs to be a permutation of 1..n-1.
    #
    # For example, size = 3 we have p=(1,2) for the circle 0->1->2->0 and p=(2,1) for the circle 0->2->1->0
    cdef array.array cycle_perm
    cdef int *p_cycle_perm
    cdef int element, cur

    if not for_caching:
        cycle_perm = int_array(size)  # Don't bother creating a new buffer each time
        p_cycle_perm = ptr(cycle_perm)

    trivial = tuple(range(1, size))   # (1,2,...size-1)
    for necklace in itertools.permutations(trivial):
        if for_caching:
            cycle_perm = int_array(size)  # Create a new buffer for each circle
            p_cycle_perm = ptr(cycle_perm)

        # The actual necklace is [0]+necklace
        cycle_perm = int_array(size)  # Create a new buffer for each circle
        p_cycle_perm = ptr(cycle_perm)

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
        return _calc_all_circle_permutations(size, for_caching=False)

    if not size in _circle_cache:
        entries = list(_calc_all_circle_permutations(size, for_caching=True))
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

    cdef int i, converted_circle

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

            p_cycle = ptr(cycle)
            for circle_perm in _all_circle_permutations(cycle_len):
                p_circle_perm = ptr(circle_perm)
                # _all_circle_permtuations yields (1, 2, 0) and (2, 0 ,1)
                # The permutations we need to return are (1, 3, 2, 0) and (3, 0, 2, 1) - these have
                # one stationary point - 2, and a cycle of length 3.
                # To do this, we need to convert the cycle from (1,2,0) and (2,0,1) to "group space": (1,3,0) and (3,0,1)
                for i in range(cycle_len):
                    # Each iteration overwrites the previous iteration value, so there's no need to save the permutation
                    # between iterations.
                    converted_circle = p_cycle[p_circle_perm[i]] # Convert circle to group space
                    p_perm[p_cycle[i]] = converted_circle
                yield from generate(cycle_index-1)    # Apply the rest of the circles

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


def molecule_permuter(int molecule_size, groups, int cycle_size, int add_cycles_of_two):
    """
    Generates all permutations of a molecule
    :param molecule_size: Molecule size
    :param groups: Equivalency groups
    :param cycle_size: Allowed cycle size
    :param add_cycles_of_two: When true, cycles of size of two are legal
    :return: Generator for all the permutations
    """

    cdef array.array perm = int_array(molecule_size)
    cdef int *p_perm = ptr(perm)

    cdef int i
    for i in range(molecule_size):
        p_perm[i] = i

    def generate(group_idx):
        """ Goes over all the permutations of the first group, applies each
         permutation and recursively aplies permutations of the entire groups
        :param group_idx: Index of current group
        :return:
        """
        cdef array.array group_array
        cdef int *p_group
        cdef int group_len

        cdef array.array group_perm
        cdef int *p_group_perm

        cdef int molecule_space, i

        if group_idx < 0:
            yield perm
            return

        group = groups[group_idx]
        group_len = len(group)
        if group_len == 1:
            yield from generate(group_idx-1)
        else:
            group_array = array.array('i', groups[group_idx])
            p_group = ptr(group_array)

            # Example:
            # Lets say the permutation is (0, 1 ,2 ,3), and the group is (0, 1, 3)

            for group_perm in group_permuter(len(group), cycle_size, add_cycles_of_two):
                p_group_perm = ptr(group_perm)

                # group_permuter yields (1, 2, 0) and (2, 0 ,1)
                # The permutations we need to return are (1, 3, 2, 0) and (3, 0, 2, 1) - these have
                # one stationary point - 2, and a cycle of length 3.
                # To do this, we need to convert the group_perm from (1,2,0) and (2,0,1) to (1,3,0) and (3,0,1),
                # we call this 'molecule_space'
                for i in range(group_len):
                    molecule_space = p_group[p_group_perm[i]]
                    p_perm[p_group[i]] = molecule_space
                yield from generate(group_idx-1)    # Apply the rest of the circles

    yield from generate(len(groups)-1)

