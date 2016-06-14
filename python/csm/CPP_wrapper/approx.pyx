from libc.math cimport sqrt
import numpy as np
cimport numpy as np
from csm.calculations.constants import MAXDOUBLE

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
        # print("Creating DistanceMatrix for group of size ", self.group_size)
        self.mv_distances = np.ones((len(group), len(group)), order="c") * MAXDOUBLE
        self._allowed_rows = np.zeros(len(group), dtype='i')
        self._allowed_cols = np.zeros(len(group), dtype='i')
        # print("DistanceMatrix created")

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
            self.myprint()
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
            # print("Can't find next in cycle. Allowed columns:")
            # for i in range(self.group_size):
            #   if self._allowed_cols[i]:
            #        print(i, '-->', searched_row[i])
            self.myprint()
            raise ValueError("Can't find next in cycle. Constraints: %s" % str(constraints))
        return (from_val, min_i)

    def myprint(self):
        print(self.mv_distances.base)


cdef class Vector3DHolder

def estimate_perm(op_type, op_order, molecule, dir,  chainperm):
    def cycle_builder(chainperm):
        used=np.zeros(len(chainperm)) #is it already in a cycle? 0-no 1-yes
        cycles=[]
        for i in range(len(chainperm)):
            cycle=[]
            if used[i]==1:
                continue
            cycle_head=chainperm[i]
            used[cycle_head]=1
            cycle.append(cycle_head)
            next_index=chainperm[cycle_head]
            while next_index!=cycle_head:
                used[next_index]=1
                cycle.append(next_index)
                next_index=chainperm[next_index]
            print("cycle", cycle)
            cycles.append(cycle)
        print("cycles", cycles)
        return cycles


    #print("estimate_perm(dir=%s) called, type of dir is %s" % (dir, type(dir)))
    # create rotation matrix
    rotation_mat = create_rotation_matrix(1, op_type, op_order, dir)
    # run rotation matrix on atoms
    rotated = (rotation_mat @ molecule.Q.T).T
    cdef Vector3DHolder rotated_holder = Vector3DHolder(rotated)

    cdef Vector3DHolder Q_holder = Vector3DHolder(molecule.Q)

    cdef double *a
    cdef double *b

    # create permutation:
    perm = [-1] * len(molecule)

    #permutation creation is done by group:
    #print("Estimating permutation for dir ", dir)
    for i in range(len(molecule.equivalence_classes)):
        for cycle in cycle_builder(chainperm):
            group=[]
            for index in cycle:
                print("chain in cycle", molecule.chain_groups[i][index])
                group = group + molecule.chain_groups[i][index]
            chain_len=len( molecule.chain_groups[i][index])
            print("group:",group)
            print("group len", len(group), "chain len", chain_len)
            distances = DistanceMatrix(group)

            for from_index in cycle:
                from_chain=molecule.chain_groups[i][from_index]
                to_chain=molecule.chain_groups[i][cycle[from_index]]
                print("from chain", from_chain, "to chain", to_chain)
                for j in from_chain:
                    for k in to_chain:
                        a = rotated_holder.get_vector(j)
                        b = Q_holder.get_vector(k)
                        distance = array_distance(a,b)
                        distances.add(group.index(k), group.index(j), distance)
        perm = perm_builder(op_type, op_order, group, distances, perm, chainperm)

    # print("Returning estimated permutation")
    return perm

def perm_builder(op_type, op_order, group, distance_matrix, perm, chainperm):
    print("Building permutation for group of len %d" % len(group))
    group_id=np.min(group)
    left=len(group)
    while left>=op_order:
        print("Building cycle (left=%d)..." % left)
        (from_val, to_val)=distance_matrix.get_min_val()
        perm[group[from_val]]=group[to_val]
        print("%d --> %d" % (from_val, to_val))
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

            print("%d --> %d" % (from_val, to_val))
            perm[group[from_val]]=group[to_val]
            distance_matrix.remove(from_val, to_val)
            left-=1
            cycle_length+=1

    #for remaining pairs or singles, simply go through them and set them
    #TODO: this section is wrong for chains
    while True:
        print("entered problem section")
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
    print("exiting perm bulder")
    return perm


# TODO: Move this back to Python (csm_calculations), or keep it in Cython and call it from Python
def create_rotation_matrix(iOp, op_type, op_order, dir):
    is_improper = op_type != 'CN'
    is_zero_angle = op_type == 'CS'
    W = np.array([[0.0, -dir[2], dir[1]], [dir[2], 0.0, -dir[0]], [-dir[1], dir[0], 0.0]])
    rot = np.zeros((3, 3))
    angle = 0.0 if is_zero_angle else 2 * np.pi * iOp / op_order
    factor = -1 if is_improper and (iOp % 2) == 1 else 1

    # The rotation matrix is calculated similarly to the Rodrigues rotation matrix. The only
    # difference is that the matrix is also a reflection matrix when factor is -1.
    #
    # This is why we took the old C++ code instead of applying the Rodrigues formula directly.
    for s in range(3):
        for t in range(3):
            ang = np.cos(angle) if s == t else 0
            rot[s][t] = ang + ((factor - np.cos(angle)) * dir[s] * dir[t] + np.sin(angle) * W[s][t])

    return rot