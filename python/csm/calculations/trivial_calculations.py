from csm.fast import CythonPermuter
from csm.calculations.basic_calculations import CSMState, process_results, Operation
from csm.calculations.constants import MAXDOUBLE
from csm.calculations.exact_calculations import ExactCalculation, Calculation
from csm.calculations.permuters import ConstraintPermuter
from csm.molecule.molecule import Molecule, MoleculeFactory

'''
contains the trivial calculation (identity perm on a chain permutation)
and the perm count
'''

class TrivialCalculation(Calculation):
    def __init__(self, operation, molecule, use_chains=True, *args, **kwargs):
        super().__init__(operation, molecule)
        self.use_chains=use_chains
        self.calc()

    def calc(self):
        molecule=self.molecule
        if molecule.chains and self.use_chains:
            best = CSMState(molecule=molecule, op_type=self.operation.type, op_order=self.operation.order, csm=MAXDOUBLE)
            chain_permutations = []
            dummy = MoleculeFactory.dummy_molecule_from_size(len(molecule.chains), molecule.chain_equivalences)
            permuter = CythonPermuter(dummy, self.operation.order, self.operation.type, keep_structure=False, precalculate=False)
            for state in permuter.permute():
                chain_permutations.append(list(state.perm))
            for chainperm in chain_permutations:
                perm = [-1 for i in range(len(molecule))]
                for f_index in range(len(chainperm)):
                    f_chain = molecule.chains[f_index]
                    t_chain = molecule.chains[chainperm[f_index]]
                    for i in range(len(f_chain)):
                        perm[f_chain[i]] = t_chain[i]

                ec = ExactCalculation(self.operation, molecule, keep_structure=False, perm=perm)
                result=ec.result
                if result.csm < best.csm:
                    best = result

        else:
            perm = [i for i in range(len(molecule))]
            ec= ExactCalculation(self.operation, molecule, keep_structure=False, perm=perm)
            best= ec.result

        self._csm_result = process_results(best)
        return self.result


def trivial_calculation(op_type, op_order, molecule, sn_max=8, use_chains=True, *args, **kwargs):
    tc=TrivialCalculation(Operation.placeholder(op_type, op_order, sn_max), molecule, use_chains)
    return tc.result

