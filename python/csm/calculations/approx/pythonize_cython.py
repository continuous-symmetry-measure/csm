import numpy as np
from csm.calculations.constants import MAXDOUBLE
from csm.calculations.basic_calculations import create_rotation_matrix
from cython_munkres import munkres

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
        return self.mv_distances.base

    def tostr(self):
        ##print(self.mv_distances.base)
        pass




def approximate_perm_hungarian(op_type, op_order, molecule, dir, chain_perm):
    #print("Inside estimate_perm, dir=%s" % dir)
    # create rotation matrix
    rotation_mat = create_rotation_matrix(1, op_type, op_order, dir)
    # run rotation matrix on atoms
    rotated = (rotation_mat @ molecule.Q.T).T
    rotated_holder = Vector3DHolder(rotated)
    Q_holder = Vector3DHolder(molecule.Q)
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
            distances=fill_distance_matrix(len(current_atom_indices), cycle, chains_in_group, chain_perm, rotated_holder, Q_holder, atom_to_matrix_indices)

            #3. call the perm builder on the group
            perm = hungarian_perm_builder(op_type, op_order, current_atom_indices, distances, perm)
            #(4. either continue to next group or finish)
    #print(perm)
    return perm


def munkres_wrapper(A):
    x = A.shape[0]
    y = A.shape[1]
    results=[]
    res_mat=munkres(A)
    for i in range(x):
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
