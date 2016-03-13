import timeit

import array

from copy import deepcopy

from CythonPlayground.playground import CalcState, Cache, one_iter, ArrayHolder, run_iters
import numpy as np

OP_ORDER = 4
MOLECULE_SIZE = 12
NUM_ITERS = 5000000


def init_state(molecule_size, op_order):
    state = CalcState(molecule_size, op_order)
    for i in range(op_order):
        rand = np.random.permutation(molecule_size)
        state.perms.set_perm(i, rand)
    return state


def run():
    cache = Cache(size=MOLECULE_SIZE)

    print("Molecule size: %d\nOp order: %d\n\nNumber of iterations: %d" % (MOLECULE_SIZE, OP_ORDER, NUM_ITERS))
    state = init_state(MOLECULE_SIZE, OP_ORDER)
    group = list(range(OP_ORDER))
    print("Here we go")
    run_iters(state, group, cache, NUM_ITERS)

def run_arrayholder():
    for i in range(NUM_ITERS):
        holder = ArrayHolder()

if __name__=='__main__':
    timer = timeit.Timer(run)
    print(timer.timeit(number=1))


