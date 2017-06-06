from libc.math cimport sqrt
from libcpp.vector cimport vector
import numpy as np
cimport numpy as np
from csm.calculations.constants import MAXDOUBLE
from csm.calculations.basic_calculations import create_rotation_matrix
from cython_munkres import munkres

cdef class Vector3D
cdef class Matrix3D

cdef array_distance(double *a, double *b):
    return sqrt(
        (a[0]-b[0])*(a[0]-b[0])
        +(a[1]-b[1])*(a[1]-b[1])
        +(a[2]-b[2])*(a[2]-b[2]))


cdef class DistanceMatrix:
    cdef int group_size
    cdef double[:,::1] mv_distances  # Make sure this is a C-order memoryview
    cdef long[:] _allowed_rows  # _allowd_rows[i] is 1 iff the row is still available for a permutation
    cdef long[:] _allowed_cols  # Respectively.

    def __init__(self, group_size):
        self.group_size = group_size
        # ##print("Creating DistanceMatrix for group of size ", self.group_size)
        self.mv_distances = np.ones((group_size, group_size), order="c") * MAXDOUBLE
        self._allowed_rows = np.zeros(group_size, dtype='long')
        self._allowed_cols = np.zeros(group_size, dtype='long')
        # ##print("DistanceMatrix created")

    def add(self, int from_val, int to_val, double distance=MAXDOUBLE):
        self.mv_distances[from_val, to_val] = distance
        self._allowed_rows[from_val] = 1
        self._allowed_cols[to_val] = 1

    def remove(self, int from_val, int to_val):
        self._allowed_rows[from_val] = 0
        self._allowed_cols[to_val] = 0

    def get_min_val(self):
        cdef int i
        cdef int j
        cdef double min = MAXDOUBLE
        cdef int min_i = -1
        cdef int min_j = -1
        cdef double tmp

        cdef double *row_ptr
        cdef long *allowed_rows = &self._allowed_rows[0]
        cdef long *allowed_cols = &self._allowed_cols[0]

        for i in range(self.group_size):
            if self._allowed_rows[i]:
                row_ptr = &self.mv_distances[i,0]
                for j in range(self.group_size):
                    if self._allowed_cols[j]:
                        tmp = row_ptr[j]
                        if tmp < min:
                            min = tmp
                            min_i = i
                            min_j = j
        return (min_i, min_j)


    def get_next_in_cycle(self, int from_val, constraints):
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
        return self.mv_distances.base

    def tostr(self):
        ##print(self.mv_distances.base)
        pass


cdef class Vector3DHolder

def cycle_builder(chainperm):
        """
        the only thing this function does is build cycles ON THE BASIS OF CHAINPERM
        if chainperm has invalid cycles, this will return invalid cycles
        :param chainperm:
        :return: a cycle in chainperm
        """
        chainperm_copy=list(chainperm)
        ##print("chainperm_copy", chainperm_copy)

        def recursive_cycle_builder(chain_head, index, cycle):
            ##print(chain_head, index, cycle)
            cycle.append(index)
            chainperm_copy[index]=-1
            if chainperm[index]==chain_head:
                ##print("yielding cycle", cycle)
                yield cycle
                cycle=[]
            else:
                yield from recursive_cycle_builder(chain_head,chainperm[index], cycle)

        for i, val in enumerate(chainperm_copy):
            ##print("I,val", i,val)
            if val==-1:
                continue
            chainperm_copy[i]=-1
            ##print("cycle head is", i)
            yield from recursive_cycle_builder(i, i,[])
        #find and return cycles
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

cdef fill_distance_matrix(len_group, cycle, chain_group, chain_perm, Vector3DHolder rotated_holder, Vector3DHolder Q_holder, long[:] matrix_indices):
    cdef double *a
    cdef double *b
    cdef DistanceMatrix distances = DistanceMatrix(len_group)
    cdef vector[int] from_chain
    cdef vector[int] to_chain
    cdef int from_chain_index
    cdef int to_chain_index
    cdef int i
    cdef int j


    for chain_index in cycle: # Use an iterator
        # Todo: Convert from_chain and to_chain into vector[ints]
        from_chain = list_to_vector_int(chain_group[chain_index])
        to_chain = list_to_vector_int(chain_group[chain_perm[chain_index]])
        for i in range(from_chain.size()):
            from_chain_index = from_chain[i]
            for j in range(to_chain.size()):
                to_chain_index = to_chain[j]
                a = rotated_holder.get_vector(from_chain_index)
                b = Q_holder.get_vector(to_chain_index)
                distance = array_distance(a,b)
                distances.add(matrix_indices[to_chain_index], matrix_indices[from_chain_index], distance)
    return distances

def approximate_perm_classic(op_type, op_order, molecule, dir, chain_perm):
    #print("Inside estimate_perm, dir=%s" % dir)
    # create rotation matrix
    rotation_mat = create_rotation_matrix(1, op_type, op_order, dir)
    # run rotation matrix on atoms
    rotated = (rotation_mat @ molecule.Q.T).T
    cdef Vector3DHolder rotated_holder = Vector3DHolder(rotated)
    cdef Vector3DHolder Q_holder = Vector3DHolder(molecule.Q)
    cdef long[:] atom_to_matrix_indices = np.ones(len(molecule), dtype='long') * -1
    # empty permutation:
    perm = [-1] * len(molecule)

    #permutation is built by "group": equivalence class, and valid cycle within chain perm (then valid exchange w/n cycle)
    for cycle in cycle_builder(chain_perm):
        # Todo: Convert cycle into a vector[int]
        for chains_in_group in molecule.groups_with_internal_chains:
            #1. create the group of atom indices we will be building a distance matrix with
            try:
                current_atom_indices=get_atom_indices(cycle, chains_in_group)
            except KeyError: #chaingroup does not have chains belonging to current cycle
                continue #continue to next chain group
            atom_to_matrix_indices=get_atom_to_matrix_indices(current_atom_indices, atom_to_matrix_indices)

            #2. within that group, go over legal switches and add their distance to the matrix
            distances=fill_distance_matrix(len(current_atom_indices), cycle, chains_in_group, chain_perm, rotated_holder, Q_holder, atom_to_matrix_indices)

            #3. call the perm builder on the group
            perm = perm_builder(op_type, op_order, current_atom_indices, distances, perm)
            #(4. either continue to next group or finish)
    #print(perm)
    return perm


def approximate_perm_hungarian(op_type, op_order, molecule, dir, chain_perm):
    #print("Inside estimate_perm, dir=%s" % dir)
    # create rotation matrix
    rotation_mat = create_rotation_matrix(1, op_type, op_order, dir)
    # run rotation matrix on atoms
    rotated = (rotation_mat @ molecule.Q.T).T
    cdef Vector3DHolder rotated_holder = Vector3DHolder(rotated)
    cdef Vector3DHolder Q_holder = Vector3DHolder(molecule.Q)
    cdef long[:] atom_to_matrix_indices = np.ones(len(molecule), dtype='long') * -1
    # empty permutation:
    perm = [-1] * len(molecule)

    #permutation is built by "group": equivalence class, and valid cycle within chain perm (then valid exchange w/n cycle)
    for cycle in cycle_builder(chain_perm):
        # Todo: Convert cycle into a vector[int]
        for chains_in_group in molecule.groups_with_internal_chains:
            #1. create the group of atom indices we will be building a distance matrix with
            try:
                current_atom_indices=get_atom_indices(cycle, chains_in_group)
            except KeyError: #chaingroup does not have chains belonging to current cycle
                continue #continue to next chain group
            atom_to_matrix_indices=get_atom_to_matrix_indices(current_atom_indices, atom_to_matrix_indices)

            #2. within that group, go over legal switches and add their distance to the matrix
            distances=fill_distance_matrix(len(current_atom_indices), cycle, chains_in_group, chain_perm, rotated_holder, Q_holder, atom_to_matrix_indices)

            #3. call the perm builder on the group
            perm = hungarian_perm_builder(op_type, op_order, current_atom_indices, distances, perm)
            #(4. either continue to next group or finish)
    #print(perm)
    return perm

@cython.boundscheck(False)
def munkres_wrapper(np.ndarray[np.double_t,ndim=2, mode="c"] A not None):
    cdef int x = A.shape[0]
    cdef int y = A.shape[1]
    results=[]
    res_mat=munkres(A) #return a matrix of booleans with True marking optimal positions
    for i in range(x): #here we convert that matrix to a list of indices
        for j in range(y):
            if res_mat[i][j]:
                results += [(i, j)]
    return results

def hungarian_perm_builder(op_type, op_order, group, distance_matrix, perm):
    matrix=distance_matrix.get_matrix()
    #m = Munkres()
    #indexes = m.compute(matrix)
    indexes=munkres_wrapper(matrix)
    for (from_val, to_val) in indexes:
        perm[group[from_val]]=group[to_val]
    return perm

def perm_builder(op_type, op_order, group, distance_matrix, perm):
    #print("Building permutation for group of len %d" % len(group))
    left=len(group)
    while left>=op_order:
        #print("Building cycle (left=%d)..." % left)
        (from_val, to_val)=distance_matrix.get_min_val()
        perm[group[from_val]]=group[to_val]
        #print("%d --> %d" % (from_val, to_val))
        distance_matrix.remove(from_val, to_val)
        left-=1
        if from_val==to_val: #cycle length 1 completed
            continue
        #otherwise, build cycle:
        cycle_head=from_val
        cycle_length=1
        cycle_done=False

        while not cycle_done:
            constraints= set()
            constraints.add(cycle_head) #prevents too-short cycle
            if cycle_length == op_order-1:
                from_val = to_val
                to_val = cycle_head
                cycle_done = True
            else:
                constraints.add(to_val)#prevent dead-end cycle
                constraints.add(from_val)  # prevent dead-endloop
                if (op_type=='SN' or op_order==2) and cycle_length<2:
                    constraints.remove(cycle_head) #remove cycle_head from constraints (also removes from_val)

                (next_from_val, next_to_val)=distance_matrix.get_next_in_cycle(to_val, constraints)
                if next_to_val==from_val: #cycle length 2 - only possible if above if(op_type==SN... was True
                    cycle_done=True
                (from_val, to_val)=(next_from_val, next_to_val)

            #print("%d --> %d" % (from_val, to_val))
            perm[group[from_val]]=group[to_val]
            distance_matrix.remove(from_val, to_val)
            left-=1
            cycle_length+=1

    #for remaining pairs or singles, simply go through them and set them
    #TODO: this section is wrong for chains
    while True:
        (from_val, to_val) = distance_matrix.get_min_val()
        if from_val==-1:  # We're finished with this group
            break
        if perm[group[from_val]] != -1: #we've finished with the group
            break
        if op_type=='SN':
            perm[group[from_val]] = group[to_val]
            perm[group[to_val]] = group[from_val]
            distance_matrix.remove(from_val, to_val)
            distance_matrix.remove(to_val, from_val)
            #print("%d --> %d" % (from_val, to_val))
        else:
            perm[group[from_val]]=group[from_val]
            #print("%d --> %d" % (from_val, from_val))
            distance_matrix.remove(from_val, from_val)
    return perm
