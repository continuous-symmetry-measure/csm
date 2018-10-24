'''
the classes used for running the approximate algorithm
'''

import multiprocessing
import operator

import datetime
from collections import OrderedDict

import numpy as np
from csm.fast import approximate_perm_classic, munkres_wrapper

from csm.calculations.approx.dirs import ClassicDirectionChooser
from csm.input_output.formatters import format_CSM
from csm.calculations.exact_calculations import ExactCalculation
from csm.calculations.basic_calculations import create_rotation_matrix, array_distance, check_perm_cycles, \
    CalculationTimeoutError
from csm.calculations.data_classes import CSMState, Operation, CSMResult
from csm.calculations.constants import MAXDOUBLE, CSM_THRESHOLD, MINDOUBLE
from csm.molecule.molecule import MoleculeFactory
from csm.fast import CythonPermuter
from csm.input_output.formatters import csm_log as print
from csm.calculations.basic_calculations import now, run_time


class _OptionalLogger:
    '''
    Utility class for having calls to a function which, if non-existent, does nothing
    '''
    def __init__(self, log_func=None):
        self._log_func = log_func
    def _log(self, *args):
        if self._log_func:
            self._log_func(*args)

class _SingleDirectionStatistics:
    # per direction, we want to store:
    # 1. every direction passed through
    # 2. every csm passed through, and percent cycle preservation
    # 3. runtime
    def __init__(self, dir):
        self.start_dir = dir
        self.dirs = []
        self.csms = []
        self.cycle_stats = []
        self._stop_reason = ""

    def append_sub_direction(self, result):
        self.dirs.append(result.dir)
        self.csms.append(result.csm)
        self.cycle_stats.append(result.num_invalid)

    @property
    def stop_reason(self):
        return self._stop_reason

    @stop_reason.setter
    def stop_reason(self, reason):
        self._stop_reason = reason

    def start_clock(self):
        self.__start_time = datetime.datetime.now()

    def end_clock(self):
        now = datetime.datetime.now()
        time_d = now - self.__start_time
        self.run_time = time_d.total_seconds()

    def __repr__(self):
        return str({
            "dirs": self.dirs,
            "csms": self.csms
        })

    @property
    def end_dir(self):
        return self.dirs[-1]

    @property
    def start_csm(self):
        return self.csms[0]

    @property
    def end_csm(self):
        return self.csms[-1]

    @property
    def num_iterations(self):
        return len(self.dirs)

    def __lt__(self, other):
        try:  # iof other doesn't have csm, we are less than them
            that_one = other.end_csm
        except:
            return True

        try:  # if we don't have csm, other is less than us
            this_one = self.end_csm
        except:
            return False

        return self.end_csm < other.end_csm

    def to_dict(self):
        try:
            return_dict= {
            "start dir":list(self.start_dir),
            "start csm":self.start_csm,
            "stop reason": self.stop_reason,
            "end dir": list(self.end_dir),
            "end csm":self.end_csm,
            "num iterations":self.num_iterations,
            "dirs":[list(dir) for dir in self.dirs],
            "csms":self.csms,
            "cycle stats":self.cycle_stats,
            "run time":self.run_time
                }
            return return_dict
        except:
            return{
            "start dir":list(self.start_dir),
            "stop reason": "was never reached"
            }

class DirectionStatisticsContainer:
    def __init__(self, initial_directions):
        self.directions_dict = OrderedDict()
        self.directions_arr=[]
        for index, dir in enumerate(initial_directions):
            self.directions_dict[tuple(dir)] = _SingleDirectionStatistics(dir)
            self.directions_arr.append(self.directions_dict[tuple(dir)])

    def __getitem__(self, key):
        return self.directions_dict[tuple(key)]

    def __setitem__(self, key, value):
        self.directions_dict[tuple(key)] = value

    def __iter__(self):
        return self.directions_dict.__iter__()

    def __str__(self):
        return str(self.directions_dict)

    def to_dict(self):
        return [{"dir":dir, "stats":self.directions_dict[dir].to_dict()} for dir in self.directions_dict]

class ApproxStatistics(DirectionStatisticsContainer):
    pass


class SingleDirApproximator(_OptionalLogger):
    def __init__(self, operation, molecule, perm_from_dir_builder, log_func=None, timeout=100,
                 max_iterations=50):
        self._log_func = log_func
        self._molecule = molecule
        self._op_type = operation.type
        self._op_order = operation.order
        self._operation=operation
        self.max_iterations = max_iterations
        self.perm_from_dir_builder = perm_from_dir_builder(operation, molecule, log_func, timeout)
        self._chain_permutations = self.perm_from_dir_builder.get_chain_perms()

    def _create_perm_from_dir(self, dir, chainperm):
        return self.perm_from_dir_builder.create_perm_from_dir(dir, chainperm)

    def calculate(self, dir):
        statistics = _SingleDirectionStatistics(dir)
        statistics.start_clock()
        best = CSMState(molecule=self._molecule, op_type=self._op_type, op_order=self._op_order, csm=MAXDOUBLE)
        self.least_invalid = CSMState(molecule=self._molecule, op_type=self._op_type, op_order=self._op_order,
                                      csm=MAXDOUBLE, num_invalid=MAXDOUBLE)
        self._log("Calculating for initial direction: ", dir)
        for chainperm in self._chain_permutations:
            if len(self._chain_permutations) > 1:
                self._log("\tCalculating for chain permutation ", chainperm)
            best_for_chain_perm = old_results = CSMState(molecule=self._molecule, op_type=self._op_type,
                                                         op_order=self._op_order,
                                                         csm=MAXDOUBLE, dir=dir)
            i = 0
            while True:
                i += 1
                self._log("\t\titeration", i, ":")

                try:
                    perm = self._create_perm_from_dir(old_results.dir, chainperm)
                    interim_results = ExactCalculation.exact_calculation_for_approx(self._operation,
                                                                                    self._molecule, perm=perm)
                except CalculationTimeoutError as e:
                    self._log("\t\titeration ", i, " Timed out after ", str(e.timeout_delta), " seconds")
                    break

                if interim_results.num_invalid < self.least_invalid.num_invalid or \
                        (interim_results.num_invalid == self.least_invalid.num_invalid
                                     and interim_results.csm < self.least_invalid.csm):
                    self.least_invalid = interim_results
                statistics.append_sub_direction(interim_results)
                # self._log("\t\t\tfound a permutation using dir", old_results.dir, "...")
                if i > 1:
                    self._log("\t\t\tthere are",
                              len(perm) - np.sum(np.array(perm) == np.array(old_results.perm)),
                              "differences between new permutation and previous permutation")
                self._log("\t\t\tthe csm for this permutation is:", str(round(interim_results.csm, 8)))
                self._log("\t\t\tthe new direction from this permutation is", interim_results.dir)
                self._log("\t\t\tthe distance between the new direction and the previous direction is:",
                          str(round(np.linalg.norm(interim_results.dir - old_results.dir), 8)))

                if interim_results.csm < best_for_chain_perm.csm:
                    best_for_chain_perm = interim_results

                # Various stop conditions for the loop, listed as multiple if statements so that the code is clearer
                if i >= self.max_iterations:
                    statistics.stop_reason = "Max iterations"
                    self._log("\t\tStopping after %d iterations" % i)
                    break
                # if i > 1 and math.fabs(old_results.csm - interim_results.csm) / math.fabs(old_results.csm) > 0.01:
                #    self._log("\t\tStopping due to CSM ratio")
                if best_for_chain_perm.csm < CSM_THRESHOLD:
                    statistics.stop_reason = "CSM below threshold"
                    self._log("\t\tStopping because the best CSM is good enough")
                    break
                if abs(np.linalg.norm(interim_results.dir - old_results.dir)) <= 0:
                    statistics.stop_reason = "No change in direction"
                    self._log("\t\tStopping because the direction has not changed")
                    break
                if interim_results.csm >= old_results.csm:  # We found a worse CSM
                    statistics.stop_reason = "No improvement in CSM"
                    self._log("\t\tStopping because CSM did not improve (worse or equal)")
                    break

                old_results = interim_results

            if best_for_chain_perm.csm < best.csm:
                best = best_for_chain_perm
                if best_for_chain_perm.csm < CSM_THRESHOLD:
                    break

        statistics.end_clock()
        return best, statistics

class _PermFromDirBuilder(_OptionalLogger):
    def __init__(self, operation, molecule, log_func, timeout):
        self.start_time=datetime.datetime.now()
        self._log_func = log_func
        self._molecule = molecule
        self.operation=operation
        self._op_type = operation.type
        self._op_order = operation.order
        self.timeout = timeout
        self._precalculate()

    def _precalculate(self):
        '''should be overridden by inheriting classes that need chain permutations'''
        self._chain_permutations=[[0]]

    def get_chain_perms(self):
        return self._chain_permutations

    def create_perm_from_dir(self, dir, chainperm):
        raise NotImplementedError


class _ChainPermsPermBuilder(_PermFromDirBuilder):
    def _calc_chain_permutations(self):
        chain_permutations = []
        dummy = MoleculeFactory.dummy_molecule_from_size(len(self._molecule.chains), self._molecule.chain_equivalences)
        permuter = CythonPermuter(dummy, self._op_order, self._op_type, precalculate=False)
        for state in permuter.permute():
            chain_permutations.append([i for i in state.perm])
        return chain_permutations

    def _precalculate(self):
        self._chain_permutations = self._calc_chain_permutations()


class _GreedyPermBuilder(_ChainPermsPermBuilder):
    '''
    This uses the Cython implementation of the classic (greedy) approximate algorithm.
    It is not optimized for molecules with many chain permutations.
    '''

    def create_perm_from_dir(self, dir, chainperm):
        return approximate_perm_classic(self._op_type, self._op_order, self._molecule, dir, chainperm)


class _HungarianPermBuilder(_ChainPermsPermBuilder):
    '''
    This uses the Hungarian (munkres) algorithm for optimization of cost matrix.
         It is not optimized for molecules with many chain permutations.
    '''

    def create_perm_from_dir(self, dir, chainperm):
        perm = self.approximate_perm_hungarian(self._op_type, self._op_order, self._molecule, dir, chainperm)
        perm = self.cookie_dough(perm)
        return perm

    class DistanceMatrix:

        def __init__(self, group_size):
            self.group_size = group_size
            # ##print("Creating DistanceMatrix for group of size ", self.group_size)
            self.mv_distances = np.ones((group_size, group_size), order="c") * MAXDOUBLE
            self._allowed_rows = np.zeros(group_size, dtype='long')
            self._allowed_cols = np.zeros(group_size, dtype='long')
            # ##print("DistanceMatrix created")

        def add(self, from_val, to_val, distance=MAXDOUBLE):
            self.mv_distances[from_val, to_val] = distance
            self._allowed_rows[from_val] = 1
            self._allowed_cols[to_val] = 1

        def remove(self, from_val, to_val):
            self._allowed_rows[from_val] = 0
            self._allowed_cols[to_val] = 0

        def get_min_val(self):
            min = MAXDOUBLE
            min_i = -1
            min_j = -1
            allowed_rows = self._allowed_rows[0]
            allowed_cols = self._allowed_cols[0]

            for i in range(self.group_size):
                if self._allowed_rows[i]:
                    row_ptr = self.mv_distances[i, 0]
                    for j in range(self.group_size):
                        if self._allowed_cols[j]:
                            tmp = row_ptr[j]
                            if tmp < min:
                                min = tmp
                                min_i = i
                                min_j = j
            return (min_i, min_j)

        def get_next_in_cycle(self, from_val, constraints):
            if not self._allowed_rows[from_val]:
                self.tostr()
                raise ValueError("get_next_in_cycle called with an unavailable row %d" % from_val)

            searched_row = self.mv_distances[from_val]
            min = MAXDOUBLE
            min_i = -1
            for i in range(self.group_size):
                if self._allowed_cols[i] and not i in constraints:
                    tmp = searched_row[i]
                    if tmp < min:
                        min, min_i = tmp, i

            if min_i == -1:
                # ##print("Can't find next in cycle. Allowed columns:")
                # for i in range(self.group_size):
                #   if self._allowed_cols[i]:
                #        ##print(i, '-->', searched_row[i])
                self.tostr()
                raise ValueError("Can't find next in cycle. Constraints: %s" % str(constraints))
            return (from_val, min_i)

        def get_matrix(self):
            return self.mv_distances

        def tostr(self):
            ##print(self.mv_distances.base)
            pass

    def cycle_builder(self, chainperm):
        """
        the only thing this function does is build cycles ON THE BASIS OF CHAINPERM
        if chainperm has invalid cycles, this will return invalid cycles
        :param chainperm:
        :return: a cycle in chainperm
        """
        chainperm_copy = list(chainperm)

        ##print("chainperm_copy", chainperm_copy)

        def recursive_cycle_builder(chain_head, index, cycle):
            ##print(chain_head, index, cycle)
            cycle.append(index)
            chainperm_copy[index] = -1
            if chainperm[index] == chain_head:
                ##print("yielding cycle", cycle)
                yield cycle
                cycle = []
            else:
                yield from recursive_cycle_builder(chain_head, chainperm[index], cycle)

        for i, val in enumerate(chainperm_copy):
            ##print("I,val", i,val)
            if val == -1:
                continue
            chainperm_copy[i] = -1
            ##print("cycle head is", i)
            yield from recursive_cycle_builder(i, i, [])
        # find and return cycles
        pass

    def get_atom_indices(self, cycle, chain_group):
        indices = []
        for chain_index in cycle:
            indices += chain_group[chain_index]
        return indices

    def get_atom_to_matrix_indices(self, atom_indices, atom_to_matrix_indices):
        i = 0
        for index in atom_indices:
            atom_to_matrix_indices[index] = i
            i += 1
        return atom_to_matrix_indices

    def fill_distance_matrix(self, len_group, cycle, chain_group, chain_perm, rotated_holder, Q_holder, matrix_indices):
        distances = self.DistanceMatrix(len_group)

        for chain_index in cycle:  # Use an iterator
            # Todo: Convert from_chain and to_chain into vector[ints]
            from_chain = chain_group[chain_index]
            to_chain = chain_group[chain_perm[chain_index]]
            for i in range(len(from_chain)):
                from_chain_index = from_chain[i]
                for j in range(len(to_chain)):
                    to_chain_index = to_chain[j]
                    a = rotated_holder[from_chain_index]
                    b = Q_holder[to_chain_index]
                    distance = array_distance(a, b)
                    distances.add(matrix_indices[to_chain_index], matrix_indices[from_chain_index], distance)
        return distances

    def hungarian_perm_builder(self, op_type, op_order, group, distance_matrix, perm):
        matrix = distance_matrix.get_matrix()
        # m = Munkres()
        # indexes = m.compute(matrix)
        indexes = munkres_wrapper(matrix)
        for (from_val, to_val) in indexes:
            perm[group[from_val]] = group[to_val]
        return perm

    def approximate_perm_hungarian(self, op_type, op_order, molecule, dir, chain_perm):
        # print("Inside estimate_perm, dir=%s" % dir)
        # create rotation matrix
        rotation_mat = create_rotation_matrix(1, op_type, op_order, dir)
        # run rotation matrix on atoms
        rotated = (rotation_mat @ molecule.Q.T).T
        atom_to_matrix_indices = np.ones(len(molecule), dtype='long') * -1
        # empty permutation:
        perm = [-1] * len(molecule)

        # permutation is built by "group": equivalence class, and valid cycle within chain perm (then valid exchange w/n cycle)
        for cycle in self.cycle_builder(chain_perm):
            # Todo: Convert cycle into a vector[int]
            for chains_in_group in molecule.groups_with_internal_chains:
                # 1. create the group of atom indices we will be building a distance matrix with
                try:
                    current_atom_indices = self.get_atom_indices(cycle, chains_in_group)
                except KeyError:  # chaingroup does not have chains belonging to current cycle
                    continue  # continue to next chain group
                atom_to_matrix_indices = self.get_atom_to_matrix_indices(current_atom_indices, atom_to_matrix_indices)

                # 2. within that group, go over legal switches and add their distance to the matrix
                distances = self.fill_distance_matrix(len(current_atom_indices), cycle, chains_in_group, chain_perm,
                                                      rotated,
                                                      molecule.Q, atom_to_matrix_indices)

                # 3. call the perm builder on the group
                perm = self.hungarian_perm_builder(op_type, op_order, current_atom_indices, distances, perm)
                # (4. either continue to next group or finish)
        # print(perm)
        return perm

    def cookie_dough(self, perm):
        falsecount, num_invalid, cycle_counts, indices_in_bad_cycles = check_perm_cycles(perm, self.operation)
        return perm


class _ManyChainsPermBuilder(_PermFromDirBuilder):
    '''
    This approximator uses the Hungarian (munkwres) algorithm to choose the optimal permutation of chains, rather than
    iterating through all possible chain permutations. It is hence more efficient for molecules with many possible chain
    permutations
    '''

    def create_perm_from_dir(self, dir, chainperm="dont care"):

        chain_len = len(self._molecule.chains[0])
        for chain in self._molecule.chains:
            test_len = len(self._molecule.chains[chain])
            if test_len != chain_len:
                raise ValueError("--many-chains currently expects all chains to be of same length")

        rotation_mat = create_rotation_matrix(1, self._op_type, self._op_order, dir)
        perm = [-1] * len(self._molecule)
        # improved use chains algorithm:
        # for a given symmetric axis:
        # M= number of fragments in molecule. 0>i>M, 0>j>M
        num_fragments = (len(self._molecule.chains))
        # A1: matrix A of size MxM: initialize with zeros
        fragment_distance_matrix = np.ones((num_fragments, num_fragments), order="c") * MAXDOUBLE
        # A2: confirm that there's same number of atoms in each equivalence group between fragment and fragment j
        # if there isn't, we leave M[i,j]=0 and continue to next ij
        # (this is already confirmed in molecule setup, so we will use the equivalent chains from mol.chain_equivalences
        for equivalent_chain_group in self._molecule.chain_equivalences:
            for i, frag_i in enumerate(equivalent_chain_group):
                for j, frag_j in enumerate(equivalent_chain_group):
                    fragment_distance_matrix[i, j] = self._get_fragment_distance(frag_i, frag_j, rotation_mat)

        # Run the hungarian algorithm on Aij, and thereby find a permutation between the fragments
        indexes = munkres_wrapper(fragment_distance_matrix)
        # C: rerun A3 on the relevant ijs
        for (i, j) in indexes:
            frag_i_groups = self._molecule.chains_with_internal_groups[i]
            frag_j_groups = self._molecule.chains_with_internal_groups[j]
            # D: put together into a full permutation
            for k, group_k in enumerate(frag_i_groups):
                group_m = frag_j_groups[k]
                indexes, group_distance_matrix = self._hungarian_on_groups(group_k, group_m, rotation_mat)
                for (from_val, to_val) in indexes:
                    perm[group_k[from_val]] = group_m[to_val]
        return perm

    def _get_fragment_distance(self, frag_i, frag_j, rotation_mat):
        total_distance = 0
        frag_i_groups = self._molecule.chains_with_internal_groups[frag_i]
        frag_j_groups = self._molecule.chains_with_internal_groups[frag_j]
        # A3: e = number of equivalence groups in fragment i, of size N1... Ne
        # 0>k>e
        for k, group_k in enumerate(frag_i_groups):
            group_m = frag_j_groups[k]
            indexes, group_distance_matrix = self._hungarian_on_groups(group_k, group_m, rotation_mat)
            # and the result (=sum of matrix members on diagonal generalized that hungarian found) ????
            # we add to A[i,j]
            for (i, j) in indexes:
                total_distance += group_distance_matrix[i, j]
        return total_distance

    def _hungarian_on_groups(self, group_k, group_m, rotation_mat):
        # matrix B of size Nk x Nk
        if not group_k:
            return [], []  # the way we built chains with internal groups, groups with no members are None
        size = len(group_k)
        group_distance_matrix = np.ones((size, size), order="c") * MAXDOUBLE
        # Vi..Vnk are the coordinate vectors of fragment i in equivalence class k
        # Wi..Wnk are the coordinate vectors of fragment j in equivalence class k
        for a, index_v in enumerate(group_k):
            for b, index_w in enumerate(group_m):
                coord_v = self._molecule.Q[index_v]
                coord_w = self._molecule.Q[index_w]
                # T is the symmetric translation appropriate to the given symmetric axis
                translated_v = rotation_mat @ coord_v
                # in position a,b of matrix B we place the square of the distance between T(Va) and Wb
                # in other words matrix B is a single block from the standard approx algorithm's distance matrix,
                # specifically the block that matches fragments i, j and equivalence group k
                distance = array_distance(translated_v, coord_w)
                group_distance_matrix[a, b] = distance
        # B: We run the hungarian algorithm on matrix B,
        indexes = munkres_wrapper(group_distance_matrix)
        return indexes, group_distance_matrix

class ApproxCalculation(_OptionalLogger):
    def __init__(self, operation, molecule, approx_algorithm='greedy',
                 log_func=lambda *args: None, selective=False, num_selected=10, *args, **kwargs):
        self.operation=operation
        self._molecule = molecule

        self._log_func = log_func

        # choose the approximator class
        if approx_algorithm == 'hungarian':
            self.perm_builder = _HungarianPermBuilder
        if approx_algorithm == 'greedy':
            self.perm_builder = _GreedyPermBuilder
        if approx_algorithm == 'many-chains':
            self.perm_builder = _ManyChainsPermBuilder

        # get the directions
        direction_chooser=ClassicDirectionChooser(molecule, operation.type, operation.order)
        self._initial_directions = direction_chooser.dirs
        if len(self._initial_directions) > 1:
            self._log("There are", len(self._initial_directions),
                      "initial directions to search for the best permutation")


        self.selective = selective
        self.num_selected = num_selected
        self.statistics = ApproxStatistics(self._initial_directions)
        self._max_iterations = 30

    def handle_chirality(self):
        pass

    def calculate(self, timeout=100, *args, **kwargs):
        self.start_time = now()
        self.timeout = timeout
        overall_stats={}
        if self.operation.type=="CH": # Chirality
            # First CS
            best_op=Operation('cs')
            best_result = self._calculate(best_op)
            if best_result.csm > MINDOUBLE:
                # Try the SN's
                for op_order in range(2, self.operation.order + 1, 2):
                    op=Operation("S"+str(op_order))
                    result = self._calculate(self.operation)
                    if result.csm < best_result.csm:
                        best_result = result
                        best_op=op
                    if best_result.csm < MINDOUBLE:
                        break
            self.operation.order=best_op.order
            self.operation.type=best_op.type
        else:
            best_result=self._calculate(self.operation)
        overall_stats["runtime"]=run_time(self.start_time)
        self.result = CSMResult(best_result, self.operation, overall_stats=overall_stats,
                                    ongoing_stats={"approx":self.statistics.to_dict()})
        return self.result

    def _calculate(self, operation):
        if operation.type == 'CI' or (operation.type == 'SN' and operation.order == 2):
            dir = [1.0, 0.0, 0.0]
            if operation.type == 'SN':
                op_msg = 'S2'
            else:
                op_msg = 'CI'
            self._log("Operation %s - using just one direction: %s" % (op_msg, dir))
            best_results = self._calculate_for_directions(operation, [dir], 1)
        else:
            if self.selective:
                self._calculate_for_directions(operation, self._initial_directions, 1)
                best_dirs = []
                sorted_csms = sorted(self.statistics.directions_arr)
                for item in sorted_csms[:self.num_selected]:
                    best_dirs.append(item.start_dir)
                    self._log("Running again on the", self.num_selected, "best directions")
                best_results = self._calculate_for_directions(operation, best_dirs, self._max_iterations)

            else:
                best_results = self._calculate_for_directions(operation, self._initial_directions, self._max_iterations)

        best_result, least_invalid = best_results
        if least_invalid.num_invalid < best_result.num_invalid:
            if least_invalid.csm <= best_result.csm:
                best_result = least_invalid
            else:
                print("(A result with better preservation of integrity of cycle lengths was found")
                print("Direction: ", least_invalid.dir, " yields a CSM of", format_CSM(least_invalid.csm),
                      "\n",(1 - (least_invalid.num_invalid / len(self._molecule))) * 100,
                      "% of the molecule's atoms are in legal cycles)")
                print("--------")
        return best_result


    def _calculate_for_directions(self, operation, dirs, max_iterations):
        best = CSMState(molecule=self._molecule, op_type=operation.type, op_order=operation.order, csm=MAXDOUBLE,
                        num_invalid=MAXDOUBLE)
        least_invalid = CSMState(molecule=self._molecule,  op_type=operation.type, op_order=operation.order,
                                 csm=MAXDOUBLE, num_invalid=MAXDOUBLE)
        single_dir_approximator = SingleDirApproximator(operation, self._molecule,
                                                        self.perm_builder, self._log,
                                                        self.timeout, max_iterations=max_iterations)
        for dir in dirs:
            best_result_for_dir, statistics = single_dir_approximator.calculate(dir)
            self.statistics[dir]=statistics
            least_invalid_for_dir = single_dir_approximator.least_invalid
            if least_invalid_for_dir.num_invalid < least_invalid.num_invalid or \
                    (least_invalid_for_dir.num_invalid == least_invalid.num_invalid and least_invalid_for_dir.csm < least_invalid.csm):
                least_invalid = least_invalid_for_dir
            if best_result_for_dir.csm < best.csm:
                best = best_result_for_dir
                if best.csm < CSM_THRESHOLD:
                    break
        return best, least_invalid

