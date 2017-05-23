import numpy as np
from csm.calculations.approx.base import Approximator
from csm.fast import approximate_perm_classic,  munkres_wrapper
from csm.fast import approximate_perm_hungarian as cython_hungarian

from csm.calculations.basic_calculations import create_rotation_matrix
from csm.calculations.constants import MAXDOUBLE


class ClassicApproximator(Approximator):
    """
    This classic approximator delegates the calculation to the old Cython implementation of approximate_perm_classic
    """

    def _approximate(self, dir, chainperm):
        return approximate_perm_classic(self._op_type, self._op_order, self._molecule, dir, chainperm)

class CythonHungarianApproximator(Approximator):
    """
    This approximator uses the Cython implementation approximate_perm_hungariam
    """

    def _approximate(self, dir, chainperm):
        return approximate_perm_hungarian(self._op_type, self._op_order, self._molecule, dir, chainperm)
        #return cython_hungarian(self._op_type, self._op_order, self._molecule, dir, chainperm)


def array_distance(a,b):
    return np.sqrt(
        (a[0]-b[0])*(a[0]-b[0])
        +(a[1]-b[1])*(a[1]-b[1])
        +(a[2]-b[2])*(a[2]-b[2]))

class DistanceMatrix:

    def __init__(self, group_size):
        self.group_size = group_size
        # ##print("Creating DistanceMatrix for group of size ", self.group_size)
        self.mv_distances = np.ones((group_size, group_size), order="c") * MAXDOUBLE
        self._allowed_rows = np.zeros(group_size, dtype='long')
        self._allowed_cols = np.zeros(group_size, dtype='long')
        # ##print("DistanceMatrix created")

    def add(self,  from_val,  to_val, distance=MAXDOUBLE):
        self.mv_distances[from_val, to_val] = distance
        self._allowed_rows[from_val] = 1
        self._allowed_cols[to_val] = 1

    def remove(self,  from_val,  to_val):
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
                row_ptr = self.mv_distances[i,0]
                for j in range(self.group_size):
                    if self._allowed_cols[j]:
                        tmp = row_ptr[j]
                        if tmp < min:
                            min = tmp
                            min_i = i
                            min_j = j
        return (min_i, min_j)


    def get_next_in_cycle(self,  from_val, constraints):
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

        if min_i==-1:
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

def cycle_builder(chainperm):
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

def get_atom_indices(cycle, chain_group):
    indices=[]
    for chain_index in cycle:
        indices+=chain_group[chain_index]
    return indices

def get_atom_to_matrix_indices(atom_indices, atom_to_matrix_indices):
    i=0
    for index in atom_indices:
        atom_to_matrix_indices[index]=i
        i+=1
    return atom_to_matrix_indices

def fill_distance_matrix(len_group, cycle, chain_group, chain_perm, rotated_holder, Q_holder, matrix_indices):
    distances = DistanceMatrix(len_group)

    for chain_index in cycle: # Use an iterator
        # Todo: Convert from_chain and to_chain into vector[ints]
        from_chain = chain_group[chain_index]
        to_chain = chain_group[chain_perm[chain_index]]
        for i in range(len(from_chain)):
            from_chain_index = from_chain[i]
            for j in range(len(to_chain)):
                to_chain_index = to_chain[j]
                a = rotated_holder[from_chain_index]
                b = Q_holder[to_chain_index]
                distance = array_distance(a,b)
                distances.add(matrix_indices[to_chain_index], matrix_indices[from_chain_index], distance)
    return distances

def hungarian_perm_builder(op_type, op_order, group, distance_matrix, perm):
    matrix=distance_matrix.get_matrix()
    #m = Munkres()
    #indexes = m.compute(matrix)
    indexes=munkres_wrapper(matrix)
    for (from_val, to_val) in indexes:
        perm[group[from_val]]=group[to_val]
    return perm

def approximate_perm_hungarian(op_type, op_order, molecule, dir, chain_perm):
    #print("Inside estimate_perm, dir=%s" % dir)
    # create rotation matrix
    rotation_mat = create_rotation_matrix(1, op_type, op_order, dir)
    # run rotation matrix on atoms
    rotated = (rotation_mat @ molecule.Q.T).T
    atom_to_matrix_indices = np.ones(len(molecule), dtype='long') * -1
    # empty permutation:
    perm = [-1] * len(molecule)

    #permutation is built by "group": equivalence class, and valid cycle within chain perm (then valid exchange w/n cycle)
    for cycle in cycle_builder(chain_perm):
        # Todo: Convert cycle into a vector[int]
        for chains_in_group in molecule.group_chains:
            #1. create the group of atom indices we will be building a distance matrix with
            try:
                current_atom_indices=get_atom_indices(cycle, chains_in_group)
            except KeyError: #chaingroup does not have chains belonging to current cycle
                continue #continue to next chain group
            atom_to_matrix_indices=get_atom_to_matrix_indices(current_atom_indices, atom_to_matrix_indices)

            #2. within that group, go over legal switches and add their distance to the matrix
            distances=fill_distance_matrix(len(current_atom_indices), cycle, chains_in_group, chain_perm, rotated, molecule.Q, atom_to_matrix_indices)

            #3. call the perm builder on the group
            perm = hungarian_perm_builder(op_type, op_order, current_atom_indices, distances, perm)
            #(4. either continue to next group or finish)
    #print(perm)
    return perm