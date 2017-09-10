from csm.fast import CythonPermuter
from csm.calculations.basic_calculations import CSMState, process_results
from csm.calculations.constants import MAXDOUBLE
from csm.calculations.exact_calculations import csm_operation
from csm.molecule.molecule import Molecule, MoleculeFactory


def trivial_calculation(op_type, op_order, molecule, use_chains=True, *args, **kwargs):
    if molecule.chains and use_chains:
        best = CSMState(molecule=molecule, op_type=op_type, op_order=op_order, csm=MAXDOUBLE)
        chain_permutations = []
        dummy = MoleculeFactory.dummy_molecule_from_size(len(molecule._chains), molecule.chain_equivalences)
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

            result = csm_operation(op_type, op_order, molecule, keep_structure=False, perm=perm)
            if result.csm < best.csm:
                best = result

    else:
        perm = [i for i in range(len(molecule))]
        best = csm_operation(op_type, op_order, molecule, keep_structure=False, perm=perm)

    return process_results(best)
