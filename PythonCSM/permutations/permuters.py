import itertools
import warnings

__author__ = 'zmbq'

warnings.warn(
    'This is the inefficient Python implementation, use CPP_Wrapper.permutations for the efficient  version',
    RuntimeWarning)


def _get_cycle_structs(perm_size, cycle_sizes):
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
            return

        for cycle_size in cycle_sizes:
            if cycle_size > len_left:
                continue
            rest = left[1:]

            if cycle_size==1:
                yield from generate(cycles, rest)
                continue

            start = (left[0],)
            for rest in itertools.combinations(rest, cycle_size-1):
                cycle = start + rest
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
    trivial = tuple(range(1, size))   # (1,2,...size-1)
    for necklace in itertools.permutations(trivial):
        # The actual necklace is [0]+necklace
        cycle_perm = [0] * size
        cur = 0
        for element in necklace:
            cycle_perm[cur] = element
            cur = element
        cycle_perm[necklace[-1]] = 0  # Add the [0] that is missing from the necklace
        yield cycle_perm

_circle_cache = {}  # perm_size->all circles of size
_CACHE_LIMIT = 10


def all_circle_permutations(size):
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
    trivial = list(range(perm_size))

    def generate(perm, cycle_index):
        # Goes over all the circles of the first cycle, apply each circle and
        # recursively generates the circles of the rest of the cycles
        # perm is built gradually, each recursive call applies one cycle
        if cycle_index < 0:
            yield tuple(perm)
            return

        cycle = cycle_struct[cycle_index]
        if len(cycle) == 1:
            yield from generate(perm, cycle_index-1)
        else:
            # Example:
            # Lets say the permutation is (0, 1 ,2 ,3), and the cycle is (0, 1, 3)
            circles = all_circle_permutations(len(cycle))
            for circle_perm in circles:
                # _all_circle_permtuations yields (1, 2, 0) and (2, 0 ,1)
                # The permutations we need to return are (1, 3, 2, 0) and (3, 0, 2, 1) - these have
                # one stationary point - 2, and a cycle of length 3.
                # To do this, we need to convert the cycle from (1,2,0) and (2,0,1) to (1,3,0) and (3,0,1)
                for i in range(len(cycle)):
                    converted_circle = cycle[circle_perm[i]]
                    perm[cycle[i]] = converted_circle  # Perm is now (1, 3, 2, 0) or (3, 0, 2, 1) <--- The converted_circle applied to (0, 1, 3)
                yield from generate(perm, cycle_index-1)    # Apply the rest of the circles

    yield from generate(trivial[:], len(cycle_struct)-1)


def group_permuter(group_size, cycle_size, add_cycles_of_two):
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


def molecule_permuter(elements, groups, cycle_size, add_cycles_of_two):
    """
    Generates all permutations of a molecule
    :param molecule_size: Molecule size
    :param groups: Equivalency groups
    :param cycle_size: Allowed cycle size
    :param add_cycles_of_two: When true, cycles of size of two are legal
    :return: Generator for all the permutations
    """
    def generate(ordered_elements, groups_left):
        """ Goes over all the permutations of the first group, applies each
         permutation and recursively aplies permutations of the entire groups
        :param ordered_elements:
        :param groups_left:
        :return:
        """
        if not groups_left:
            yield tuple(ordered_elements)
            return

        group = groups_left[0]
        groups_left = groups_left[1:]
        if len(group) == 1:
            yield from generate(ordered_elements, groups_left)
        else:
            # Example:
            # Lets say the permutation is (0, 1 ,2 ,3), and the group is (0, 1, 3)
            start_elements_order = ordered_elements[:]
            for group_perm in group_permuter(len(group), cycle_size, add_cycles_of_two):
                # group_permuter yields (1, 2, 0) and (2, 0 ,1)
                # The permutations we need to return are (1, 3, 2, 0) and (3, 0, 2, 1) - these have
                # one stationary point - 2, and a cycle of length 3.
                # To do this, we need to convert the group_perm from (1,2,0) and (2,0,1) to (1,3,0) and (3,0,1)
                converted_group = [group[i] for i in group_perm]  # Converted circle is now (1,3,0) or (3,0,1)
                for i in range(len(converted_group)):
                    ordered_elements[group[i]] = start_elements_order[converted_group[i]]
                    # Perm is now (1, 3, 2, 0) or (3, 0, 2, 1) <--- The converted_circle applied to (0, 1, 3)
                yield from generate(ordered_elements, groups_left)    # Apply the rest of the circles
                ordered_elements = start_elements_order[:]

    elements_order = elements  # The starting elements order
    yield from generate(elements_order, groups)


print(list(molecule_permuter(['a', 'b', 'c', 'd', 'e'], [[0, 1, 3], [2, 4]], 2, False)))
