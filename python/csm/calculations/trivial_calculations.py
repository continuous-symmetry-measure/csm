import datetime
from csm.fast import CythonPermuter
from csm.calculations.data_classes import CSMState, CSMResult
from csm.calculations.constants import MAXDOUBLE
from csm.calculations.exact_calculations import ExactCalculation
from csm.molecule.molecule import MoleculeFactory
from csm.calculations.basic_calculations import run_time

'''
contains the trivial calculation (identity perm on a chain permutation)
'''

class TrivialCalculation:
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
        self.operation=operation
        self.molecule=molecule
        self.use_chains=use_chains
        self.statistics={}
        self.start_time=datetime.datetime.now()

    def calculate(self, *args, **kwargs):
        molecule=self.molecule
        if molecule.chains and self.use_chains:
            best = CSMState(molecule=molecule, op_type=self.operation.type, op_order=self.operation.order, csm=MAXDOUBLE)
            chain_permutations = []
            dummy = MoleculeFactory.dummy_molecule_from_size(len(molecule.chains), molecule.chain_equivalences)
            permuter = CythonPermuter(dummy, self.operation.order, self.operation.type, precalculate=False)
            for state in permuter.permute():
                chain_permutations.append(list(state.perm))
            for chainperm in chain_permutations:
                perm = [-1 for i in range(len(molecule))]
                for f_index in range(len(chainperm)):
                    f_chain = molecule.chains[f_index]
                    t_chain = molecule.chains[chainperm[f_index]]
                    for i in range(len(f_chain)):
                        perm[f_chain[i]] = t_chain[i]

                result= ExactCalculation.exact_calculation_for_approx(self.operation, molecule, perm=perm)
                if result.csm < best.csm:
                    best = result

        else:
            perm = [i for i in range(len(molecule))]
            best= ExactCalculation.exact_calculation_for_approx(self.operation, molecule, perm=perm)

        self._csm_result = CSMResult(best, self.operation, overall_stats={"runtime":run_time(self.start_time)})
        return self.result

    @property
    def result(self):
        return self._csm_result