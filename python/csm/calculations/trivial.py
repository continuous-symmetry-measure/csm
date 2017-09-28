from csm.fast import CythonPermuter
from csm.calculations.basic_calculations import CSMState, process_results
from csm.calculations.constants import MAXDOUBLE
from csm.calculations.exact_calculations import exact_calculation
from csm.calculations.permuters import ConstraintPermuter
from csm.molecule.molecule import Molecule, MoleculeFactory


def trivial_calculation(op_type, op_order, molecule, use_chains=True, *args, **kwargs):
    if molecule.chains and use_chains:
        best = CSMState(molecule=molecule, op_type=op_type, op_order=op_order, csm=MAXDOUBLE)
        chain_permutations = []
        dummy = MoleculeFactory.dummy_molecule_from_size(len(molecule.chains), molecule.chain_equivalences)
        permuter = CythonPermuter(dummy, op_order, op_type, keep_structure=False, precalculate=False)
        for state in permuter.permute():
            chain_permutations.append(list(state.perm))
        for chainperm in chain_permutations:
            perm = [-1 for i in range(len(molecule))]
            for f_index in range(len(chainperm)):
                f_chain = molecule.chains[f_index]
                t_chain = molecule.chains[chainperm[f_index]]
                for i in range(len(f_chain)):
                    perm[f_chain[i]]=t_chain[i]

            result = exact_calculation(op_type, op_order, molecule, keep_structure=False, perm=perm)
            if result.csm < best.csm:
                best = result

    else:
        perm = [i for i in range(len(molecule))]
        best = exact_calculation(op_type, op_order, molecule, keep_structure=False, perm=perm)

    return process_results(best)



def perm_count(op_type, op_order, molecule, keep_structure, print_perms=False, no_constraint=False, *args, **kwargs):
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

        def log_fact(n):
            """
            Returns the log(n!)
            """
            result = 0.0
            for i in range(1,n+1):
                result += math.log(i)
            return result

        def log_nCr(n, r):
            """
            Calculates the log of n choose r
            """
            if r > n:
                raise ValueError("r must be <= n for nCr")
            if r < 0:
                raise ValueError("r must be non-negative")
            if n <= 0:
                raise ValueError("n must be positive")

            return log_fact(n) - log_fact(n-r) - log_fact(r)

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


    if not print_perms:
        if not keep_structure:
            count=1
            groups=molecule.equivalence_classes
            for group in groups:
                count*=_len_group_permuter(len(group), op_order, op_type=='SN')
            return int(count)


    traced_state = CSMState(molecule=molecule, op_type=op_type, op_order=op_order)
    permuter = ConstraintPermuter(molecule, op_order, op_type, keep_structure, precalculate=False)
    if no_constraint:
        permuter=CythonPermuter(molecule, op_order, op_type, keep_structure, precalculate=False)

    for state in permuter.permute():
        if csm_state_tracer_func:
            traced_state=traced_state._replace(csm = '', perm = state.perm, dir = '')
            csm_state_tracer_func(traced_state)
        if permuter.count%1000000==0:
            print("counted", int(permuter.count/1000000), "million permutations thus far...")
        count=permuter.count
    return count
