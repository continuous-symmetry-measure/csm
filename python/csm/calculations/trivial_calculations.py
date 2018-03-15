from csm.fast import CythonPermuter
from csm.calculations.data_classes import CSMState, process_results, Operation, Calculation
from csm.calculations.constants import MAXDOUBLE
from csm.calculations.exact_calculations import ExactCalculation
from csm.calculations.permuters import ConstraintPermuter
from csm.molecule.molecule import Molecule, MoleculeFactory

'''
contains the trivial calculation (identity perm on a chain permutation)
and the perm count
'''

class TrivialCalculation(Calculation):
    """
    Calculates the CSM of the identity permutation of a molecule. 
    If use-chains is specified, calculates the identity permutation of every possible chain permutation, returns best.
    """
    def __init__(self, operation, molecule, use_chains=True, *args, **kwargs):
        """
        :param operation: instance of Operation class or named tuple, with fields for name and order, that describes the symmetry.
        :param molecule: instance of Molecule class on which the described symmetry calculation will be performed.
        :param use_chains: default True. When True, all possible chain permutations with an identity perm on their components are measured.
                When false, only the pure identity perm is measured.
        """
        super().__init__(operation, molecule)
        self.use_chains=use_chains
        self.statistics={}

    def calculate(self):
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
                ec.calculate()
                result=ec.result
                if result.csm < best.csm:
                    best = result

        else:
            perm = [i for i in range(len(molecule))]
            ec= ExactCalculation(self.operation, molecule, keep_structure=False, perm=perm)
            ec.calculate()
            best= ec.result

        self._csm_result = best
        return self.result