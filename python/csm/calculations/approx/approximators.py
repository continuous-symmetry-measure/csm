import operator

import numpy as np
import math
from csm.fast import approximate_perm_classic,  munkres_wrapper
from csm.fast import approximate_perm_hungarian as cython_hungarian
from csm.calculations.exact_calculations import exact_calculation
from csm.calculations.basic_calculations import create_rotation_matrix, array_distance, CSMState
from csm.calculations.constants import MAXDOUBLE
from csm.calculations.permuters import ApproxConstraintPermuter
from csm.molecule.molecule import Molecule, MoleculeFactory
from csm.fast import CythonPermuter


class Approximator:
    '''
    This is the base class all approximators inherit from.
    It implements an algorithm whereby initial directions are chosen (via _choose_initial_direction)
    And then iterated through in 'approximate' with the function '_approximate_from_initial_direction'
    All inheriting classes must implement _approximate_from_initial_direction, and may optionally implement _precalculate
    '''
    def __init__(self, op_type, op_order, molecule, dir_chooser, log_func=lambda *args: None, timeout=300):
        self._op_type = op_type
        self._op_order = op_order
        self._molecule = molecule
        self._initial_directions=dir_chooser.dirs
        self._log=log_func
        self._log("There are", len(self._initial_directions), "initial directions to search for the best permutation")
        self._timeout = timeout

    def approximate(self):
        # the basic steps of direction-based approximation are as follows:
        # 0. precalculation of any variables that only need to be calculated once
        self._precalculate()
        best = CSMState(molecule=self._molecule, op_type=self._op_type, op_order=self._op_order, csm=MAXDOUBLE)
        # 1. choose an initial direction
        for dir in self._initial_directions:
            #calculate on the basis of that permutation as detailed in the function
            result= self._approximate_from_initial_dir(dir)
            print("Csm is:", result.csm)
            # 5. repeat from 1, using a different starting direction (assuming more than one)
            if result.csm<best.csm:
                best=result
        # 6. return the best result
        return best

    def _precalculate(self):
        pass

    def _approximate_from_initial_dir(self, dir):
        # 2. using the direction, create a permutation
        # 3. using the permutation and it's results, choose another direction
        # 4. repeat 2-3 until an endpoint
        raise NotImplementedError


class OldApproximator(Approximator):
    '''
    This uses the Cython implementation of the classic (greedy) approximate algorithm.
    It is not optimized for molecules with many chain permutations.
    '''
    def _calc_chain_permutations(self):
        chain_permutations = []
        dummy = MoleculeFactory.dummy_molecule_from_size(len(self._molecule.chains), self._molecule.chain_equivalences)
        permuter = CythonPermuter(dummy, self._op_order, self._op_type, keep_structure=False, precalculate=False)
        for state in permuter.permute():
            chain_permutations.append([i for i in state.perm])

        return chain_permutations

    def _precalculate(self):
        self._chain_permutations = self._calc_chain_permutations()

    def _for_inversion(self, best):
        # if inversion:
        # not necessary to calculate dir, use geometrical center of structure
        dir = [1.0, 0.0, 0.0]
        if self._op_type == 'SN':
            op_msg = 'S2'
        else:
            op_msg = 'CI'
        self._log("Operation %s - using just one direction: %s" % (op_msg, dir))

        for chainperm in self._chain_permutations:
            self._log("Calculating for chain permutation ", chainperm)
            perm = self._build_perm(dir, chainperm)
            best_for_chain_perm = exact_calculation(self._op_type, self._op_order, self._molecule, keep_structure=False,
                                                perm=perm)
            if best_for_chain_perm.csm < best.csm:
                best = best_for_chain_perm

        return best

    def _approximate_from_initial_dir(self, dir):
        best = CSMState(molecule=self._molecule, op_type=self._op_type, op_order=self._op_order, csm=MAXDOUBLE)

        if self._op_type == 'CI' or (self._op_type == 'SN' and self._op_order == 2):
            return self._for_inversion(best)

        self._log("Calculating for initial direction: ", dir)
        for chainperm in self._chain_permutations:
                self._log("\tCalculating for chain permutation ", chainperm)
                # find permutation for this direction of the symmetry axis
                perm = self._build_perm(dir, chainperm)
                # solve using this perm until it converges:
                old_results = CSMState(molecule=self._molecule, op_type=self._op_type, op_order=self._op_order,
                                       csm=MAXDOUBLE)
                best_for_chain_perm = interim_results = exact_calculation(self._op_type, self._op_order, self._molecule,
                                                                      keep_structure=False, perm=perm)
                self._log("\t\tfound initial permutation")
                self._log("\t\tfirst pass yielded dir", interim_results.dir,
                            "and CSM " + str(round(interim_results.csm, 5)))
                # print(perm)

                if best_for_chain_perm.csm < best.csm:
                    best = best_for_chain_perm

                # iterations:
                i = 0
                max_iterations = 50
                while (i < max_iterations
                       and math.fabs(old_results.csm - interim_results.csm) > 0.0001
                       and interim_results.csm < old_results.csm
                       and interim_results.csm > 0.0001):
                    old_results = interim_results
                    i += 1
                    perm = self._build_perm(interim_results.dir, chainperm)
                    interim_results = exact_calculation(self._op_type, self._op_order, self._molecule, keep_structure=False,
                                                    perm=perm)

                    self._log("\t\titeration", i, ":")
                    self._log("\t\t\tfound a permutation using dir", old_results.dir, "...")
                    self._log("\t\t\tthere are",
                              len(perm) - np.sum(np.array(perm) == np.array(old_results.perm)),
                                "differences between new permutation and previous permutation")
                    self._log("\t\t\tusing new permutation, found new direction", interim_results.dir)
                    self._log("\t\t\tthe distance between the new direction and the previous direction is:",
                              str(round(np.linalg.norm(interim_results.dir - old_results.dir), 8)))
                    self._log("\t\t\tthe csm found is:", str(round(interim_results.csm, 8)))
                    # print(perm)

                    if interim_results.csm < best_for_chain_perm.csm:
                        diff = best_for_chain_perm.csm - interim_results.csm
                        best_for_chain_perm = interim_results
                        if best_for_chain_perm.csm < best.csm:
                            best = best_for_chain_perm

        return best

    def _build_perm(self, dir, chainperm):
        return approximate_perm_classic(self._op_type, self._op_order, self._molecule, dir, chainperm)


class HungarianApproximator(OldApproximator):
    '''
    This uses the Hungarian (munkres) algorithm for optimization of cost matrix.
         It is not optimized for molecules with many chain permutations.
    '''
    def _build_perm(self, dir, chainperm):
        return self.approximate_perm_hungarian(self._op_type, self._op_order, self._molecule, dir, chainperm)

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


class ManyChainsApproximator(Approximator):
    '''
    This approximator uses the Hungarian (munkwres) algorithm to choose the optimal permutation of chains, rather than
    iterating through all possible chain permutations. It is hence more efficient for molecules with many possible chain 
    permutations
    '''
    def _approximate_from_initial_dir(self, dir):
        best = CSMState(molecule=self._molecule, op_type=self._op_type, op_order=self._op_order, csm=MAXDOUBLE)
        old_results = CSMState(molecule=self._molecule, op_type=self._op_type, op_order=self._op_order,
                               csm=MAXDOUBLE)

        perm = self._build_perm(dir)
        interim_results = exact_calculation(self._op_type, self._op_order, self._molecule,
                      keep_structure=False, perm=perm)

        # iterations:
        i = 0
        max_iterations = 50
        while (i < max_iterations and
                   (math.fabs(old_results.csm - interim_results.csm) / math.fabs(
                       old_results.csm) > 0.01
                    and interim_results.csm < old_results.csm)
                    and interim_results.csm > 0.0001):
            old_results = interim_results
            i += 1
            perm = self._build_perm(interim_results.dir)
            interim_results = exact_calculation(self._op_type, self._op_order, self._molecule, keep_structure=False,
                                            perm=perm)

            self._log("\t\titeration", i, ":")
            self._log("\t\t\tfound a permutation using dir", old_results.dir, "...")
            self._log("\t\t\tthere are",
                      len(perm) - np.sum(np.array(perm) == np.array(old_results.perm)),
                        "differences between new permutation and previous permutation")
            self._log("\t\t\tusing new permutation, found new direction", interim_results.dir)
            self._log("\t\t\tthe distance between the new direction and the previous direction is:",
                      str(round(np.linalg.norm(interim_results.dir - old_results.dir), 8)))
            self._log("\t\t\tthe csm found is:", str(round(interim_results.csm, 8)))

            if interim_results.csm < best.csm:
                best =interim_results

        return best

    def _build_perm(self, dir):
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
                    fragment_distance_matrix[i,j]=self._get_fragment_distance(frag_i, frag_j, rotation_mat)

        # Run the hungarian algorithm on Aij, and thereby find a permutation between the fragments
        indexes = munkres_wrapper(fragment_distance_matrix)
        # C: rerun A3 on the relevant ijs
        for (i,j) in indexes:
            frag_i_groups=self._molecule.chains_with_internal_groups[i]
            frag_j_groups=self._molecule.chains_with_internal_groups[j]
        # D: put together into a full permutation
            for k, group_k in enumerate(frag_i_groups):
                group_m=frag_j_groups[k]
                indexes, group_distance_matrix= self._hungarian_on_groups(group_k, group_m, rotation_mat   )
                for (from_val, to_val) in indexes:
                    perm[group_k[from_val]] = group_m[to_val]
        return perm


    def _get_fragment_distance(self, frag_i, frag_j, rotation_mat):
        total_distance=0
        frag_i_groups = self._molecule.chains_with_internal_groups[frag_i]
        frag_j_groups = self._molecule.chains_with_internal_groups[frag_j]
        # A3: e = number of equivalence groups in fragment i, of size N1... Ne
        # 0>k>e
        for k, group_k in enumerate(frag_i_groups):
            group_m=frag_j_groups[k]
            indexes, group_distance_matrix= self._hungarian_on_groups(group_k, group_m, rotation_mat)
            # and the result (=sum of matrix members on diagonal generalized that hungarian found) ????
            # we add to A[i,j]
            for (i, j) in indexes:
                total_distance+=group_distance_matrix[i,j]
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


class StructuredApproximator(OldApproximator):
    def _approximate_from_initial_dir(self, dir):
        best = CSMState(molecule=self._molecule, op_type=self._op_type, op_order=self._op_order, csm=MAXDOUBLE)
        old_results = CSMState(molecule=self._molecule, op_type=self._op_type, op_order=self._op_order,
                               csm=MAXDOUBLE, dir=dir)

        try:
            perm = self._build_perm(dir)
            best = interim_results = exact_calculation(self._op_type, self._op_order, self._molecule,
                          keep_structure=False, perm=perm)
            # iterations:
            i = 0
            max_iterations = 50
            while (i < max_iterations and
                       (math.fabs(old_results.csm - interim_results.csm) / math.fabs(
                           old_results.csm) > 0.01
                        and interim_results.csm < old_results.csm)
                        and interim_results.csm > 0.0001
                   and abs(np.linalg.norm(interim_results.dir - old_results.dir))>0):
                old_results = interim_results
                i += 1
                perm = self._build_perm(interim_results.dir)
                interim_results = exact_calculation(self._op_type, self._op_order, self._molecule, keep_structure=False,
                                                perm=perm)

                self._log("\t\titeration", i, ":")
                self._log("\t\t\tfound a permutation using dir", old_results.dir, "...")
                self._log("\t\t\tthere are",
                          len(perm) - np.sum(np.array(perm) == np.array(old_results.perm)),
                            "differences between new permutation and previous permutation")
                self._log("\t\t\tusing new permutation, found new direction", interim_results.dir)
                self._log("\t\t\tthe distance between the new direction and the previous direction is:",
                          str(round(np.linalg.norm(interim_results.dir - old_results.dir), 8)))
                self._log("\t\t\tthe csm found is:", str(round(interim_results.csm, 8)))

                if interim_results.csm < best.csm:
                    best =interim_results

                return best

        except TimeoutError:
            return best



    def _build_perm(self, dir):
        return self.build_perm_and_state(self._op_type, self._op_order, self._molecule, dir)

    def build_perm_and_state(self, op_type, op_order, molecule, dir):
        print("dir is:", dir)
        rotation_mat = create_rotation_matrix(1, op_type, op_order, dir)
        rotated = (rotation_mat @ molecule.Q.T).T
        distances_dict={}
        distances_list=[]
        for index_a, a in enumerate(molecule.Q):
            distances_dict[index_a]={}
            for index_b, b in enumerate(rotated):
                if index_b in molecule.atoms[index_a].equivalency:
                    distance = array_distance(a, b)
                else:
                    distance=MAXDOUBLE
                distances_dict[index_a][index_b]=distance
                distances_list.append(((index_a, index_b), distance))
        distances_list.sort(key=operator.itemgetter(1))
        permuter=ApproxConstraintPermuter(self._molecule, self._op_order, self._op_type, distances_dict, distances_list,
                                          timeout=self._timeout)

        state=permuter.permute().__next__()
        perm=state.perm
        return perm



