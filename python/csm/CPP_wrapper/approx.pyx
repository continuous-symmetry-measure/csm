from libc.math cimport sqrt
import numpy as np
cimport numpy as np
from csm.calculations.constants import MAXDOUBLE
from csm.calculations.basic_calculations import create_rotation_matrix

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
    cdef int[:] _allowed_rows  # _allowd_rows[i] is 1 iff the row is still available for a permutation
    cdef int[:] _allowed_cols  # Respectively.

    def __init__(self, group):
        self.group_size = len(group)
        # #print("Creating DistanceMatrix for group of size ", self.group_size)
        self.mv_distances = np.ones((len(group), len(group)), order="c") * MAXDOUBLE
        self._allowed_rows = np.zeros(len(group), dtype='i')
        self._allowed_cols = np.zeros(len(group), dtype='i')
        # #print("DistanceMatrix created")

    def add(self, int from_val, int to_val, double distance=MAXDOUBLE):
        self.mv_distances[from_val, to_val] = distance
        self._allowed_rows[from_val] = 1
        self._allowed_cols[to_val] = 1

    def sort(self):
        pass#self.list_distances = newlist = sorted(self.list_distances, key=lambda x: x.distance)

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
        cdef int *allowed_rows = &self._allowed_rows[0]
        cdef int *allowed_cols = &self._allowed_cols[0]

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

    def get_min_val_old(self):
        a = self.mv_distances
        argmin=np.argmin(a)
        from_val = int(argmin/ len(a[0]))
        to_val= int(argmin % len(a[0]))
        return (from_val, to_val)

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
            # #print("Can't find next in cycle. Allowed columns:")
            # for i in range(self.group_size):
            #   if self._allowed_cols[i]:
            #        #print(i, '-->', searched_row[i])
            self.tostr()
            raise ValueError("Can't find next in cycle. Constraints: %s" % str(constraints))
        return (from_val, min_i)

    def tostr(self):
        #print(self.mv_distances.base)
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
        #print("chainperm_copy", chainperm_copy)

        def recursive_cycle_builder(chain_head, index, cycle):
            #print(chain_head, index, cycle)
            cycle.append(index)
            chainperm_copy[index]=-1
            if chainperm[index]==chain_head:
                #print("yielding cycle", cycle)
                yield cycle
                cycle=[]
            else:
                yield from recursive_cycle_builder(chain_head,chainperm[index], cycle)

        for i, val in enumerate(chainperm_copy):
            #print("I,val", i,val)
            if val==-1:
                continue
            chainperm_copy[i]=-1
            #print("cycle head is", i)
            yield from recursive_cycle_builder(i, i,[])
        #find and return cycles
        pass

def estimate_perm(op_type, op_order, molecule, dir,  chainperm, use_chains):
    print("Inside estimate_perm, dir=%s" % dir)

    # create rotation matrix
    rotation_mat = create_rotation_matrix(1, op_type, op_order, dir)
    # run rotation matrix on atoms
    rotated = (rotation_mat @ molecule.Q.T).T
    #print("Chainperm", chainperm)
    cdef Vector3DHolder rotated_holder = Vector3DHolder(rotated)
    cdef Vector3DHolder Q_holder = Vector3DHolder(molecule.Q)
    cdef double *a
    cdef double *b

    # empty permutation:
    perm = [-1] * len(molecule)

    #permutation is built by "group": equivalence class, and valid cycle within chain perm (then valid exchange w/n cycle)
    for cycle in cycle_builder(chainperm):
        print("cycle", cycle)
        for chain_group in molecule.chain_groups:
            #print("chain group", chain_group)
            #1. find the "group" we will be building a distance matrix with
            group=[]
            for chain_index in cycle:
                group+=chain_group[chain_index]
            print("group of length %d" % len(group))
            distances = DistanceMatrix(group)

            #2. within that group, go over legal switches and add their distance to the matrix
            for chain_index in cycle:
                from_chain=chain_group[chain_index]
                to_chain=chain_group[chainperm[chain_index]]
                for j in from_chain:
                     for k in to_chain:
                         a = rotated_holder.get_vector(j)
                         b = Q_holder.get_vector(k)
                         distance = array_distance(a,b)
                         distances.add(group.index(k), group.index(j), distance)
            print("Distance matrix complete")
        #3. call the perm builder on the group
            perm = perm_builder(op_type, op_order, group, distances, perm, chainperm)
            print("Built part of the perm")
            #(4. either continue to next group or finish)

    print(perm)
    return perm

def perm_builder(op_type, op_order, group, distance_matrix, perm, chainperm):
    #print("Building permutation for group of len %d" % len(group))
    #group_id=np.min(group)
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
        else:
            perm[group[from_val]]=group[from_val]
            distance_matrix.remove(from_val, from_val)
    return perm
