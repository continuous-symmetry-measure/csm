import os
import pytest
from csm.molecule.molecule import Molecule

inbal_folder=r'C:\Users\devora.CHELEM\Sources\csm\test_cases\inbal\reading fragments'
#test for the following formats:
#pdb without chains
#pdb with chains
#pdb with hetatm
#pdb with use-sequence
os.path.join(inbal_folder, 'water-6.pdb')
#pdb with hetatm no flag
os.path.join(inbal_folder, 'water-6.pdb')
#mol with fragments
os.path.join(inbal_folder,'water-6.mol')
#mol with fragments but no fragment flag
os.path.join(inbal_folder,'water-6.mol')
#mol without fragments
r'C:\Users\devora.CHELEM\Sources\csm\test_cases\original\test16\tetralin.mol'
#csm with chains
#csm without chains
#xyz
r'C:\Users\devora.CHELEM\Sources\csm\test_cases\original\test1\AgCu10p1.xyz'
#cif


def test_pdb_with_hetatm():
    file=os.path.join(inbal_folder, 'water-6.pdb')
    Molecule.from_file(file)



#test for the following calculation types:
#approx with use-chains
#approx without use-chains
#trivial with use-chains
#trivial without use-chains
#exact with keep-structure
#normalizations with use-chains, 0-6