import pytest
from os import path
from csm.molecule.molecule import Molecule, MoleculeFactory
from csm.input_output.arguments import get_split_arguments
from csm.main.csm_run import run as csmrun
from conftest import test_folder

def xtest_pdb_with_chains_hetatm_and_nonequal_equivalence_classes():
    file=r'D:\UserData\devora\Sources\csm\test_cases\old_test_cases\inbal\proteins\2xql.pdb'
    m=MoleculeFactory.from_file(file, format='pdb', use_chains=True)


    args= "--approx --use-chains --remove-hy --many-chains"

def xtest_bonds_exist():
    '''
    read connectivity didn't call create bond set and as a result the check for existent bonds failed.
    '''
    pass

#File "C:\Users\devora.CHELEM\Sources\csm\python\csm\calculations\approx\approximators.py", line 347, in _approximate
#indexes, group_distance_matrix= self._hungarian_on_groups(group_k, group_m, rotation_mat   )
#File "C:\Users\devora.CHELEM\Sources\csm\python\csm\calculations\approx\approximators.py", line 386, in _hungarian_on_group_group_distance_matrix[a, b] = distance
#IndexError: index 121 is out of bounds for axis 1 with size 121

def xtest_use_sequence_shouldnt_auto_use_babel():
    path1=path.join(test_folder, '2hyn-noh001.pdb')
    path2=path.join(test_folder, '2hyn-noh050.pdb')
    args=get_split_arguments('c5', path1, "out.txt", "--approx", "--use-sequence")
    m = MoleculeFactory.from_file(**args)


def test_index_error_with_pdb_use_sequence():
    #previously was getting an index error on molecule creation
    mol_file=path.join(test_folder, "model-endmdl-withids.pdb")
    args = get_split_arguments('c5', mol_file, "out.txt", "--approx", "--use-sequence")
    m = MoleculeFactory.from_file(**args)
