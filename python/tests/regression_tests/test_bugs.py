import pytest
from csm.molecule.molecule import Molecule
from csm.main.csm_run import run as csmrun


def test_pdb_with_chains_hetatm_and_nonequal_equivalence_classes():
    file=r'D:\UserData\devora\Sources\csm\test_cases\old_test_cases\inbal\proteins\2xql.pdb'
    m=Molecule.from_file(file, format='pdb', use_chains=True)


    args= "--approx --use-chains --remove-hy --many-chains"


#File "C:\Users\devora.CHELEM\Sources\csm\python\csm\calculations\approx\approximators.py", line 347, in _approximate
#indexes, group_distance_matrix= self._hungarian_on_groups(group_k, group_m, rotation_mat   )
#File "C:\Users\devora.CHELEM\Sources\csm\python\csm\calculations\approx\approximators.py", line 386, in _hungarian_on_group_group_distance_matrix[a, b] = distance
#IndexError: index 121 is out of bounds for axis 1 with size 121
