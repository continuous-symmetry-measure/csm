import numpy as np
import math as m
from datetime import datetime




def now():
    return datetime.now()

def run_time(start_time):
    now = datetime.now()
    time_d = now - start_time
    return time_d.total_seconds()

class CalculationTimeoutError(TimeoutError):
    def __init__(self, timeout_delta, *args, **kwargs):
        super().__init__("Calculation timed out after "+str(timeout_delta)+" seconds", *args, **kwargs)
        self.timeout_delta=timeout_delta

def check_timeout(local_start, local_timeout):
    from csm.calculations.constants import global_start_time, global_time_out
    runtime=run_time(local_start)
    if runtime>local_timeout:
        raise CalculationTimeoutError(runtime)
    g_runtime=run_time(global_start_time)
    if g_runtime>global_time_out:
        raise CalculationTimeoutError(g_runtime)


def cart2sph(x,y,z, normalize=True):
    def do_normalize(x, y, z):
        norm= x*x + y*y + z*z
        return x/norm, y/norm, z/norm
    #https://stackoverflow.com/questions/4116658/faster-numpy-cartesian-to-spherical-coordinate-conversion
    if normalize:
        x, y, z = do_normalize(x, y, z)
    XsqPlusYsq = x**2 + y**2
    r = m.sqrt(XsqPlusYsq + z**2)               # r
    elev = m.atan2(z,m.sqrt(XsqPlusYsq))     # theta
    az = m.atan2(y,x)                           # phi
    return r, elev, az




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






def check_perm_cycles(perm, operation):
    '''
    This function checks the cycles in a given permutation according to the provided operation order and type.
    It counts the legal and illegal cycles in the permutation
    :param perm:
    :param op_order: 
    :param op_type: 
    :return: the number of illegal cycles, the number of molecules in illegal cycles, a dictionary of cycle lengths-
    with key =length cycle, val= mnumber of cycles of that length, and an array of the indices in bad cycles
    '''
    op_order=operation.order
    op_type=operation.type
    checked=[False]*len(perm)
    num_molecules_in_bad_cycles=0
    num_good_cycles=0
    num_bad_cycles=0

    cycle_counts={}

    indices_in_bad_cycles=[]


    for i, index in enumerate(perm):
        if checked[i]:
            continue
        checked[i]=True
        cycle=[index]
        while not checked[index]:
            checked[index] = True
            index = perm[index]
            cycle.append(index)

        cycle_len=len(cycle)

        if cycle_len in cycle_counts:
            cycle_counts[cycle_len]+=1
        else:
            cycle_counts[cycle_len]=1

        if cycle_len == 1 or cycle_len == op_order or (cycle_len == 2 and op_type == 'SN'):
            num_good_cycles += 1
        else:
            num_molecules_in_bad_cycles += cycle_len
            num_bad_cycles+=1
            indices_in_bad_cycles+=cycle


    return num_bad_cycles, num_molecules_in_bad_cycles, cycle_counts, indices_in_bad_cycles






def check_perm_equivalence(mol, perm):
    for origin, destination in enumerate(perm):
        if destination not in mol.atoms[origin].equivalency:
            return False
    return True

def check_perm_structure_preservation(mol, perm):
    '''
    checks what percent of a permutation does not break bonds in the molecule--
    so a return value of 1 would be a perfectly preserved bond structure, of 0.5 would mean half of the molecules bonds are broken
    :param mol: the molecule being permuted
    :param perm: the permutation whose structure preservation is being measured
    :return: percent of bonds that are preserved
    '''
    if len(mol.bondset)==0:
        raise ValueError("Molecule does not have any bond information")

    broken=0
    for origin, destination in enumerate(perm):
        for adjacent in mol.atoms[origin].adjacent:
            if (destination, perm[adjacent]) not in mol.bondset:
                broken+=1

    percent_structure= (len(mol.bondset)- broken )/len(mol.bondset)

    return percent_structure


def array_distance(a, b):
    return np.sqrt(
        (a[0] - b[0]) * (a[0] - b[0])
        + (a[1] - b[1]) * (a[1] - b[1])
        + (a[2] - b[2]) * (a[2] - b[2]))