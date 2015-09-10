import itertools
import math
from permutations.utils import log_nCr, log_fact

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
            start_perm = perm[:]
            circles = _all_circle_permutations(len(cycle))
            for circle_perm in circles:
                # _all_circle_permtuations yields (1, 2, 0) and (2, 0 ,1)
                # The permutations we need to return are (1, 3, 2, 0) and (3, 0, 2, 1) - these have
                # one stationary point - 2, and a cycle of length 3.
                # To do this, we need to convert the cycle from (1,2,0) and (2,0,1) to (1,3,0) and (3,0,1)
                converted_circle = [cycle[i] for i in circle_perm]  # Converted circle is now (1,3,0) or (3,0,1)
                for i in range(len(converted_circle)):
                    perm[cycle[i]] = converted_circle[i]  # Perm is now (1, 3, 2, 0) or (3, 0, 2, 1) <--- The converted_circle applied to (0, 1, 3)
                yield from generate(perm, cycle_index-1)    # Apply the rest of the circles
                perm = start_perm[:]

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

def _len_group_permuter(group_size, cycle_size, add_cycles_of_two):
    """
    Returns the length of the group permuter.
    :param group_Size: Group size
    :param cycle_size: Cycle size
    :param add_cycles_of_two: When true, cycles of size of two are also legal
    :return: Number of permutations that will be returned by the group_permuter
    This is done by enumerating the different *kinds* of cycle_structs, and calculating the number of permutations
    for each kind.
    A *kind* of cycle struct is: (n1 cycles of size c1, n2 cycles of size c2, ...). Where c1..cn are either 1, cycle_size
    or 2 (if add_cycles_of_two is True), and the sum of ni*ci adds to group_size
    """

    def generate_kinds():
        """
        Generates all cycle_struct kinds with the given group_size and cycle_lengths.
        A kind is a mapping from cycle_size to the number of cycles
        """
        add_two = add_cycles_of_two and cycle_size!=2  # Do we add cycles of two?
        for num_cycle_size in range(0, group_size // cycle_size + 1):
            if add_two:
                left = group_size - cycle_size * num_cycle_size  # amount left
                for num_two in range(0, left // 2 + 1):
                    yield {cycle_size: num_cycle_size, 2: num_two}
            else:
                yield {cycle_size: num_cycle_size}

    def log_structs_of_kind(kind):
        """
        Returns the log of the number of cycle_structs of a given kind.
        :param kind: The kind of cycle_struct
        :return: The number of cycle_structs of kind.
        The calculation is done by choosing the cycles.

        For example, if there's only one cycle of size 2, there are group_size choose 2 different structs of this kind.
        If there are two cycles of size 2, there are group_size choose 2 * (group_size-2) choose 2 / 2
        The reasoning here is - first choose the first cycle, then choose the second cycle. Divide by two because the
        order of the two cycles doesn't matter.

        For 3 cycles of size 4 we get:
        group_size choose 4 * (group_size-4) choose 4 * (group_size-8) choose 4 / 3!
        We choose the first group of 4, then the second, then the third, but divide by 3! because the order of groups
        doesn't matter.

        If the group size is 12 and we allow 2 cycles of size 3 and 3 cycles of size 2, we get
        (12C3 * 9C3 / 2) * (6C2 * 4C2 * 2C2) / 2!
        We first choose the groups of size 3, then the groups of size 2.

        The above expressions can be further simplified, but we don't do it - this verison is more readable and performs
        very well.
        """
        log_result = 0
        size_left = group_size
        for cycle_size, cycle_count in kind.items():
            for i in range(cycle_count):
                log_result += log_nCr(size_left, cycle_size)
                size_left -= cycle_size
            log_result -= log_fact(cycle_count)
        return log_result

    def log_perms_in_kind(kind):
        """
        Returns the log of the number of permutations of a cycle_struct of a specific kind
        The number of permutations for a cycle of size n is (n-1)! so we just go and add that (add, because
        we're returning the log
        """
        log_result = 0
        for cycle_size, cycle_count in kind.items():
            log_result += cycle_count * log_fact(cycle_size-1)
        return log_result

    count = 0
    for kind in generate_kinds():
        log_num_structs = log_structs_of_kind(kind)
        log_num_perms = log_perms_in_kind(kind)
        count += math.exp(log_num_structs + log_num_perms)

    return count


def molecule_permuter(molecule_size, groups, cycle_size, add_cycles_of_two):
    """
    Generates all permutations of a molecule
    :param molecule_size: Molecule size
    :param groups: Equivalency groups
    :param cycle_size: Allowed cycle size
    :param add_cycles_of_two: When true, cycles of size of two are legal
    :return: Generator for all the permutations
    """
    def generate(perm, groups_left):
        """ Goes over all the permutations of the first group, applies each
         permutation and recursively aplies permutations of the entire groups
        :param perm:
        :param groups_left:
        :return:
        """
        if not groups_left:
            yield tuple(perm)
            return

        group = groups_left[0]
        groups_left = groups_left[1:]
        if len(group) == 1:
            yield from generate(perm, groups_left)
        else:
            # Example:
            # Lets say the permutation is (0, 1 ,2 ,3), and the group is (0, 1, 3)
            start_perm = perm[:]
            for group_perm in group_permuter(len(group), cycle_size, add_cycles_of_two):
                # group_permuter yields (1, 2, 0) and (2, 0 ,1)
                # The permutations we need to return are (1, 3, 2, 0) and (3, 0, 2, 1) - these have
                # one stationary point - 2, and a cycle of length 3.
                # To do this, we need to convert the group_perm from (1,2,0) and (2,0,1) to (1,3,0) and (3,0,1)
                converted_group = [group[i] for i in group_perm]  # Converted circle is now (1,3,0) or (3,0,1)
                for i in range(len(converted_group)):
                    perm[group[i]] = converted_group[i]  # Perm is now (1, 3, 2, 0) or (3, 0, 2, 1) <--- The converted_circle applied to (0, 1, 3)
                yield from generate(perm, groups_left)    # Apply the rest of the circles
                perm = start_perm[:]

    perm = list(range(molecule_size))  # The trivial permutation
    yield from generate(perm, groups)



def len_molecule_permuter(molecule, op_order, op_type):
    result = 1
    for group in molecule.equivalence_classes:
        result *= _len_group_permuter(len(group), op_order, op_type == 'SN')
    return result

