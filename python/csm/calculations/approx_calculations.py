import numpy as np
import math
import logging
from collections import namedtuple
from csm.calculations.constants import MINDOUBLE, MAXDOUBLE
from csm.fast import CythonPermuter, SinglePermPermuter, TruePermChecker, PQPermChecker, CythonPIP
from csm.fast import external_get_eigens as cppeigen
from csm.calculations.csm_calculations import csm_operation, CSMState, create_rotation_matrix, process_results
from csm.molecule.molecule import Molecule

logger = logging.getLogger("csm")


def trivial_calculation(op_type, op_order, molecule, use_chains=True, *args, **kwargs):
    if molecule.chains and use_chains:
        best = CSMState(molecule=molecule, op_type=op_type, op_order=op_order, csm=MAXDOUBLE)
        chainkeys=list(molecule.chains.keys())
        chain_permutations = []
        dummy = Molecule.dummy_molecule(len(molecule._chains))
        permuter = CythonPermuter(dummy, op_order, op_type, TruePermChecker, perm_class=CythonPIP)
        for state in permuter.permute():
            chain_permutations.append(list(state.perm))
        for chainperm in chain_permutations:
            perm = [-1 for i in range(len(molecule))]
            for f_index in range(len(chainperm)):
                f_chain_key=chainkeys[f_index]
                t_chain_key=chainkeys[chainperm[f_index]]
                f_chain=molecule.chains[f_chain_key]
                t_chain = molecule.chains[t_chain_key]
                for i in range(len(f_chain)):
                    perm[f_chain[i]]=t_chain[i]

            result = csm_operation(op_type, op_order, molecule, SinglePermPermuter,
                           TruePermChecker, perm)
            if result.csm < best.csm:
                best = result

    else:
        perm = [i for i in range(len(molecule))]
        best = csm_operation(op_type, op_order, molecule, SinglePermPermuter,
                               TruePermChecker, perm)

    return process_results(best)

def is_legal(perm, molecule):
    legal=True

    switch=True
    for key in molecule.chains:
        if 0 in molecule.chains[key] and perm[0] in molecule.chains[key]:
            switch=False
        break


    switched=True
    for i in range(len(perm)):
        if i in molecule.chains[key] and perm[i] in molecule.chains[key]:
            switched=False
        else:
            switched=True
        if switched != switch:
            test=perm[i]
            first= molecule.atoms[i].chain
            second= molecule.atoms[test].chain
            legal=False
            break

    return legal



def approx_calculation(op_type, op_order, molecule, detect_outliers=False,use_chains=False, *args, **kwargs):
    results= find_best_perm(op_type, op_order, molecule, detect_outliers, use_chains)
    #print(is_legal(results.perm, molecule))
    return process_results(results)


def find_best_perm(op_type, op_order, molecule, detect_outliers, use_chains):
    if op_type == 'CI' or (op_type == 'SN' and op_order == 2):
        # if inversion:
        # not necessary to calculate dir, use geometrical center of structure
        dir = [1.0, 0.0, 0.0]
        perm = estimate_perm(op_type, op_order, molecule, dir)
        best = csm_operation(op_type, op_order, molecule, SinglePermPermuter, TruePermChecker, perm)

    else:
        best = CSMState(molecule=molecule, op_type=op_type, op_order=op_order, csm=MAXDOUBLE)

        chain_permutations = []
        if molecule.chains and use_chains:
            chainkeys = list(molecule.chains.keys())
            dummy = Molecule.dummy_molecule(len(molecule.chains))
            permuter = CythonPermuter(dummy, op_order, op_type, TruePermChecker, perm_class=CythonPIP)
            for state in permuter.permute():
                chain_permutations.append([i for i in state.perm])
        else:
            chain_permutations.append([])

        dirs = find_symmetry_directions(molecule, detect_outliers, op_type)
        for dir in dirs:
            for chainperm in chain_permutations:
                if chainperm==chain_permutations[0]:
                    continue
                # find permutation for this direction of the symmetry axis
                perm = estimate_perm(op_type, op_order, molecule, dir, chainperm)

                # solve using this perm until it converges:
                best_for_this_dir = interim_results = csm_operation(op_type, op_order, molecule, SinglePermPermuter,
                                                                    TruePermChecker, perm)
                old_results = CSMState(molecule=molecule, op_type=op_type, op_order=op_order, csm=MAXDOUBLE)
                i = 0
                max_iterations = 50
                while (i < max_iterations and math.fabs(
                            old_results.csm - interim_results.csm) > 0.01 and interim_results.csm > 0.0001):
                    old_results = interim_results
                    i += 1
                    perm = estimate_perm(op_type, op_order, molecule, interim_results.dir, chainperm)
                    interim_results = csm_operation(op_type, op_order, molecule, SinglePermPermuter, TruePermChecker, perm)
                    if interim_results.csm < best_for_this_dir.csm:
                        best_for_this_dir = interim_results
                #print("attempt for dir" + str(dir) + ": best csm is:" + str(best_for_this_dir.csm) + " after " + str(i) + " iterations")


                if best_for_this_dir.csm < best.csm:
                    best = best_for_this_dir

    return best


min_group_for_outliers = 10


def find_symmetry_directions(molecule, detect_outliers, op_type):
    # get average position of each equivalence group:
    group_averages = []
    for group in molecule.equivalence_classes:
        sum = [0, 0, 0]
        for index in group:
            sum += molecule.Q[index]
        average = sum / len(group)
        group_averages.append(average)
    group_averages = np.array(group_averages)

    dirs = dir_fit(group_averages)
    if detect_outliers and len(molecule.equivalence_classes) > min_group_for_outliers:
        dirs = dirs_without_outliers(dirs, group_averages, op_type)
    dirs = dirs_orthogonal(dirs)
    return dirs


def normalize_dir(dir):
    dir = np.array(dir)  # lists also get sent to this function, but can't be /='d
    norm = math.sqrt(math.fabs(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]))
    dir /= norm
    return dir

def dir_fit(positions):
    # get vector of the average position
    sum = np.einsum('ij->j', positions)
    average = sum / len(positions)

    mat = np.zeros((3, 3), dtype=np.float64, order="c")
    for pos in positions:
        for j in range(3):
            for k in range(3):
                mat[j][k] += (pos[j] - average[j]) * (pos[k] - average[k])

    # computer eigenvalues and eigenvectors
    # lambdas - list of 3 eigenvalues of matrix (function demands them, but we won't be using them here)
    # dirs- list of 3 eigenvectors of matrix
    dirs = np.zeros((3, 3), dtype=np.float64, order="c")
    lambdas = np.zeros((3), dtype=np.float64, order="c")
    cppeigen(mat, dirs, lambdas)

    # normalize result:
    dirs = np.array([normalize_dir(dir) for dir in dirs])

    return dirs

def dirs_without_outliers(dirs, positions, op_type):
    def compute_distance_from_line(group_avg_point, test_dir_end):
        def magnitude(point1, point2):
            vec = point2 - point1
            return math.sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2])

        def get_U(Point, LineStart, LineEnd, LineMag):
            return (((Point[0] - LineStart[0]) * (LineEnd[0] - LineStart[0])) +
                    ((Point[1] - LineStart[1]) * (LineEnd[1] - LineStart[1])) +
                    ((Point[2] - LineStart[2]) * (LineEnd[2] - LineStart[2]))) / (LineMag * LineMag)

        def get_intersect(LineStart, LineEnd, U):
            Intersection = [0, 0, 0]
            Intersection[0] = LineStart[0] + U * (LineEnd[0] - LineStart[0]);
            Intersection[1] = LineStart[1] + U * (LineEnd[1] - LineStart[1]);
            Intersection[2] = LineStart[2] + U * (LineEnd[2] - LineStart[2]);
            return Intersection

        linestart = [0, 0, 0]
        linemag = magnitude(test_dir_end, linestart)
        U = get_U(group_avg_point, linestart, test_dir_end, linemag)
        intersection = get_intersect(linestart, test_dir_end, U)
        return magnitude(group_avg_point, intersection)
    print("======================detecting outliers============================")
    more_dirs = []
    for dir in dirs:
        # 1.Find the distance of each point from the line/plane
        dists=[]
        for i in range(len(positions)):
            pos = positions[i]
            if op_type == 'CS':
                dist = math.abs(dir[0] * pos[0] + dir[1] * pos[1] + dir[2] * pos[2])
            else:
                dist = compute_distance_from_line(pos, dir)
            dists.append(dist)

        # 2. Find median of the distances m
        median = np.median(np.array(dists))
        # 3. for each distance di, if di > 2*m, remove it as an outlier
        with_outliers_removed = list()
        for i in range(len(positions)):
            pos = positions[i]
            dist = dists[i]
            if not (dist / median > 2 or dist / median > 2):
                with_outliers_removed.append(pos)
        # 4. recompute dirs
        outliers_dirs = dir_fit(with_outliers_removed)
        more_dirs += list(outliers_dirs)
    return np.array(more_dirs)

def dirs_orthogonal(dirs):
    added_dirs = list(dirs)

    for dir in dirs:
        if math.fabs(dir[0]) < MINDOUBLE:
            dir1 = [1.0, 0.0, 0.0]
            dir2 = [0.0, -dir[2], dir[1]]
        elif math.fabs(dir[1]) < MINDOUBLE:
            dir1 = [-dir[2], 0.0, dir[0]]
            dir2 = [0.0, 1.0, 0.0]
        else:
            dir1 = [-dir[1], dir[0], 0.0]
            dir2 = [0.0, -dir[2], dir[1]]
        # normalize dir1:
        dir1 = normalize_dir(dir1)
        # remove projection of dir1 from dir2
        scal = dir1[0] * dir2[0] + dir1[1] * dir2[1] + dir1[2] * dir2[2]
        dir2[0] -= scal * dir1[0]
        dir2[1] -= scal * dir1[1]
        dir2[2] -= scal * dir1[2]
        # normalize dir2
        dir2 = normalize_dir(dir2)
        added_dirs.append(dir1)
        added_dirs.append(dir2)

    return np.array(added_dirs)


class distance_record:
    def __init__(self, row, col, dist):
        self.row = row
        self.distance = dist
        self.col = col

    def __str__(self):
        return "distance" + str(self.row) + "," + str(self.col) + ": " + str(self.distance)

Distance= namedtuple('Distance', 'from_val to_val distance')

class DistanceMatrix:
    def __init__(self, group):
        self.np_distances= np.ones((len(group), len(group))) * MAXDOUBLE
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
        a=self.np_distances
        from_val = int(np.argmin(a) / len(a[0]))
        to_val= int(np.argmin(a) % len(a[0]))
        return (from_val, to_val)

    def get_next_in_cycle(self, from_val, constraints):
        searched_row=self.np_distances[from_val]
        indices= np.argsort(searched_row)
        i=0
        if constraints:
            while indices[i] in constraints:
                i+=1
        to_val=indices[i]
        return (from_val, to_val)



def estimate_perm(op_type, op_order, molecule, dir, chainperm=[]):
    # this is a naive implementation, which works well with the naive c++ build_perm
    # step two: first estimation of the symmetry measure: first permutation
    # apply the symmetry operation once, using symm_element from step one

    # create rotation matrix
    rotation_mat = create_rotation_matrix(1, op_type, op_order, dir)
    # run rotation matrix on atoms
    rotated = (rotation_mat @ molecule.Q.T).T

    # create permutation:
    perm = [-1] * len(molecule)

     #permutation creation is done by group:
    for group in molecule.equivalence_classes:
        # measure all the distances between all the points in X to all the points in Y
        distances = DistanceMatrix(group)#list()
        for i in range(len(group)):
            for j in range(len(group)):
                 #legal chainperm permutation
                distance = np.linalg.norm(rotated[group[i]] - molecule.Q[group[j]])
                if chainperm:
                    number_index_of_chain_from = molecule.chainkeys.index(molecule.atoms[group[i]].chain)
                    number_index_of_chain_to = molecule.chainkeys.index(molecule.atoms[group[j]].chain)
                    if chainperm[number_index_of_chain_from] != number_index_of_chain_to:
                        distance = MAXDOUBLE
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
        #print(group_id, left, group[from_val], group[to_val], sep=", ")
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
            #print(group_id, left, group[from_val], group[to_val], sep=", ")
            left-=1
            cycle_length+=1
        pass

    #for remaining pairs or singles, simply go through them and set them
    while True:
        (from_val, to_val) = distance_matrix.get_min_val()
        if perm[group[from_val]] != -1: #we've finished with the group
            break
        if op_type=='SN':
            perm[group[from_val]] = group[to_val]
            perm[group[to_val]] = group[from_val]
            distance_matrix.remove(from_val, to_val)
            distance_matrix.remove(to_val, from_val)
            #print(group_id, left, group[from_val], group[to_val], sep=", ")
            #left-=1
            #if from_val!=to_val:
            #    print(group_id, left, group[to_val], group[from_val], sep=", ")
            #    left-=1
        else:
            perm[group[from_val]]=group[from_val]
            distance_matrix.remove(from_val, from_val)
            #left -= 1
            #print(group_id, left, group[from_val], group[from_val], sep=", ")
    return perm
