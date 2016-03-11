import timeit
from copy import deepcopy

from CythonPlayground.playground import CalcState, Factory, Cache, one_iter
import numpy as np

OP_ORDER = 6
MOLECULE_SIZE = 15
NUM_ITERS = 100000


def init_state(molecule_size, op_order, factory):
    state = CalcState(molecule_size, op_order, factory)
    for i in range(op_order):
        rand = np.random.permutation(molecule_size)
        state.perms[i, :] =rand
    return state


def run():
    factory = Factory()
    cache = Cache(size=MOLECULE_SIZE, factory=factory)

    print("Molecule size: %d\nOp order: %d\n\nNumber of iterations: %d" % (MOLECULE_SIZE, OP_ORDER, NUM_ITERS))
    state = init_state(MOLECULE_SIZE, OP_ORDER, factory)
    group = list(range(OP_ORDER))
    for i in range(NUM_ITERS):
        old_state = deepcopy(state)
        one_iter(state, group, cache)
        state = old_state

if __name__=='__main__':
    timer = timeit.Timer(run)
    print(timer.timeit(number=1))


