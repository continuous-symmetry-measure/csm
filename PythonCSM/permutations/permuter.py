import itertools

__author__ = 'zmbq'

def _get_cycle_structs(perm_size, cycle_sizes):
    """
    Generates a list of cycles in a permutation. The cycles cover the entire permutation, and are only of sizes in cycle_sizes
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
                break
            start = [left[0]]
            for rest in itertools.combinations(left[1:], cycle_size-1):
                cycle = start + list(rest)
                new_left = [item for item in left if item not in cycle]
                yield from generate(cycles+[cycle], new_left)

    if not 1 in cycle_sizes:
        raise ValueError("1 must be a valid cycle size")
    yield from generate([], list(range(perm_size)))

def _all_circles(elements):
    """ Generates all the circles (permutations that are full cycles) of 'elements'

    If elements are [A,B,C], returns (B,C,A) and (C,A,B)
    """

    def all_circle_indices(size):
        """ Returns the indices of the cycle, that will later be applied to elements """
        # To compute a full cycle of length n, we take a permutation p of size n-1 and
        # create the cycle like so: 0 goes to p[0], p[0] to p[1], p[1] to p[2] and so on, until p[n-1] goes back to 0
        # For this to work, p needs to be a permutation of 1..n-1.
        #
        # For example, size = 3 we have p=(1,2) for the circle 0->1->2->0 and p=(2,1) for the circle 0->2->1->0
        trivial = tuple(range(1, size))
        for necklace in itertools.permutations(trivial):
            # The actual necklace is [0]+necklace
            full_cycle = [0] * size
            cur = 0
            for element in necklace:
                full_cycle[cur] = cur = element
            full_cycle[necklace[-1]] = 0  # Add the [0] that is missing from the necklace
            yield full_cycle

    # yield elements  # Trivial circle
    for cycle_indices in all_circle_indices(len(elements)):
        circle = tuple(elements[i] for i in cycle_indices)
        yield circle

def _all_perms_from_cycle_struct(perm_size, cycle_struct):
    """
    Generates all permutations with a specific cycle structure
    :param perm_size: Permutation size
    :param cycle_struct: Cycle structure (list of cycles, each cycle is itself a list of its indices)
    :return: Generator that generates all the permutations
    """
    trivial = list(range(perm_size))

    def generate(perm, cycle_struct, cycles_left):
        # Goes over all the circles of the first cycle, apply each circle and
        # recursively generates the circles of the rest of the cycles
        if not cycles_left:
            yield tuple(perm)
            return

        cycle = cycles_left[0]
        cycles_left = cycles_left[1:]
        if len(cycle) > 1:
            start_perm = perm[:]
            for circle in _all_circles(cycle):
                # Apply each circle
                for i in range(len(circle)):
                    perm[cycle[i]] = circle[i]
                yield from generate(perm, cycle_struct, cycles_left) # Apply the rest of the circles
                perm = start_perm[:]
        else:
            yield from generate(perm, cycle_struct, cycles_left)

    yield from generate(trivial[:], cycle_struct, cycle_struct)


def permuter(perm_size, cycle_size, add_cycles_of_twp):
    """
    Generates all permutations of size and cycle sizes.
    :param perm_size: Permutation size
    :param cycle_size: Cycle size
    :param add_cycles_of_twp: When true, cycles of size of two are legal
    :return: Generator for all the permutations
    """
    cycle_lengths = {1, cycle_size}
    if add_cycles_of_twp:
        cycle_lengths.add(2)
    cycle_lengths = sorted(cycle_lengths)

    for cycle_struct in _get_cycle_structs(perm_size, cycle_lengths):
        # Loop over all cycle structures, which are of the allowed lengths covering the entire permutation
        yield from _all_perms_from_cycle_struct(perm_size, cycle_struct) # Return all permutations for each cycle structure
