from csm.calculations.approx.base import Approximator
from csm.fast import approximate_perm_classic, approximate_perm_hungarian

class ClassicApproximator(Approximator):
    """
    This classic approximator delegates the calculation to the old Cython implementation of approximate_perm_classic
    """

    def _approximate(self, dir, chainperm):
        return approximate_perm_classic(self._op_type, self._op_order, self._molecule, dir, chainperm)

class CythonHungarianApproximator(Approximator):
    """
    This approximator uses the Cython implementation approximate_perm_hungariam
    """

    def _approximate(self, dir, chainperm):
        return approximate_perm_hungarian(self._op_type, self._op_order, self._molecule, dir, chainperm)
