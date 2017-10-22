from csm.fast import CythonPermuter
from csm.calculations.basic_calculations import CSMState, process_results
from csm.calculations.constants import MAXDOUBLE
from csm.calculations.exact_calculations import exact_calculation
from csm.calculations.permuters import ConstraintPermuter
from csm.molecule.molecule import Molecule, MoleculeFactory

'''
contains the trivial calculation (identity perm on a chain permutation)
and the perm count
'''

def trivial_calculation(op_type, op_order, molecule, use_chains=True, *args, **kwargs):
    """
    Calculates the CSM of the identity permutation of a molecule. 
    If use-chains is specified, calculates the identity permutation of every possible chain permutation, returns best
    :param op_type: type of symmetry (CS, CN, CH, CI, SN)
    :param op_order: order of symmetry (2, 3, 4...)
    :param molecule: instance of Molecule class whose symmetry is being measured
    :param use_chains: 
    :return: 
    """
    if molecule.chains and use_chains:
        best = CSMState(molecule=molecule, op_type=op_type, op_order=op_order, csm=MAXDOUBLE)
        chain_permutations = []
        dummy = MoleculeFactory.dummy_molecule_from_size(len(molecule.chains), molecule.chain_equivalences)
        permuter = CythonPermuter(dummy, op_order, op_type, keep_structure=False, precalculate=False)
        for state in permuter.permute():
            chain_permutations.append(list(state.perm))
        for chainperm in chain_permutations:
            perm = [-1 for i in range(len(molecule))]
            for f_index in range(len(chainperm)):
                f_chain = molecule.chains[f_index]
                t_chain = molecule.chains[chainperm[f_index]]
                for i in range(len(f_chain)):
                    perm[f_chain[i]]=t_chain[i]

            result = exact_calculation(op_type, op_order, molecule, keep_structure=False, perm=perm)
            if result.csm < best.csm:
                best = result

    else:
        perm = [i for i in range(len(molecule))]
        best = exact_calculation(op_type, op_order, molecule, keep_structure=False, perm=perm)

    return process_results(best)
