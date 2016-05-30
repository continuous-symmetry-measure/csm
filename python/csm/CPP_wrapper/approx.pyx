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
    cdef np_distances
    cdef double[:,:] mv_distances
    cdef double[:] _np_delete

    def __init__(self, group):
        self.group_size = len(group)
        self.np_distances = np.ones((len(group), len(group)), order="c") * MAXDOUBLE
        self.mv_distances = self.np_distances
        self._np_delete = np.copy(self.np_distances[0])
        #self.list_distances=[]

    def add(self, from_val, to_val, distance=MAXDOUBLE):
        self.np_distances[from_val][to_val]=distance
        #dist=Distance(from_val=from_val, to_val=to_val, distance=distance)
        #self.list_distances.append(dist)

    def sort(self):
        pass#self.list_distances = newlist = sorted(self.list_distances, key=lambda x: x.distance)

    def remove(self, from_val, to_val):
        self.np_distances[:,to_val]=self._np_delete
        self.np_distances[from_val,:]=self._np_delete

    def get_min_val(self):
        cdef int i
        cdef int j
        cdef double min = MAXDOUBLE
        cdef int min_i
        cdef int min_j
        cdef double tmp

        for i in range(self.group_size):
            for j in range(self.group_size):
                tmp = self.mv_distances[i][j]
                if tmp < min:
                    min = tmp
                    min_i = i
                    min_j = j
        return (min_i, min_j)
        # a=self.np_distances
        # argmin=np.argmin(a)
        # from_val = int(argmin/ len(a[0]))
        # to_val= int(argmin % len(a[0]))
        # return (from_val, to_val)

    def get_next_in_cycle(self, from_val, constraints):
        searched_row=self.np_distances[from_val]
        indices= np.argsort(searched_row)
        i=0
        if constraints:
            while indices[i] in constraints:
                i+=1
        to_val=indices[i]
        return (from_val, to_val)


cdef class Vector3DHolder

def estimate_perm(op_type, op_order, molecule, dir, chainperm=[]):
    # this is a naive implementation, which works well with the naive c++ build_perm
    # step two: first estimation of the symmetry measure: first permutation
    # apply the symmetry operation once, using symm_element from step one

    # create rotation matrix
    rotation_mat = create_rotation_matrix(1, op_type, op_order, dir)
    # run rotation matrix on atoms
    rotated = (rotation_mat @ molecule.Q.T).T
    cdef Vector3DHolder rotated_holder = Vector3DHolder(rotated)

    #cdef double *t = rotated_holder.get_vector(1000)
    #print("Rotated holder: %f, %f, %f" % (t[0], t[1], t[2]))
    #print("Rotated: %f, %f, %f" % (rotated[1000][0], rotated[1000][1], rotated[1000][2]))
    cdef Vector3DHolder Q_holder = Vector3DHolder(molecule.Q)

    cdef double *a
    cdef double *b

    # create permutation:
    perm = [-1] * len(molecule)

     #permutation creation is done by group:
    for i in range(len(molecule.equivalence_classes)):
        group=molecule.equivalence_classes[i]
        distances = DistanceMatrix(group)
        if chainperm:
            chain_group=molecule.chain_groups[i]
            chain_indices=molecule.chain_indices[i]
            for index in chainperm:
                from_chain=chain_group[index]
                to_chain=chain_group[chainperm[index]]
                for j in from_chain:
                    for k in to_chain: #these should be the same length, but helps keep track of what j/k are
                        a = rotated_holder.get_vector(j)
                        b = Q_holder.get_vector(k)
                        distance = array_distance(a,b)
                        distances.add(chain_indices[k], chain_indices[j], distance)
        else:
            for i in range(len(group)):
                for j in range(len(group)):
                    a = rotated_holder.get_vector(i)
                    b = Q_holder.get_vector(j)
                    distance = array_distance(a,b)
                    distances.add(j,i,distance)

        perm = perm_builder(op_type, op_order, group, distances, perm, chainperm)

    return perm

def perm_builder(op_type, op_order, group, distance_matrix, perm, chainperm):
    group_id=np.min(group)
    left=len(group)
    while left>=op_order:
        (from_val, to_val)=distance_matrix.get_min_val()
        perm[group[from_val]]=group[to_val]
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
            if cycle_length== op_order-1:
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

            perm[group[from_val]]=group[to_val]
            distance_matrix.remove(from_val, to_val)
            left-=1
            cycle_length+=1
        pass

    #for remaining pairs or singles, simply go through them and set them
    #TODO: this section is wrong for chains
    while True:
        (from_val, to_val) = distance_matrix.get_min_val()
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