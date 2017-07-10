import os
import pytest
from csm.molecule.molecule import Molecule

inbal_folder=r'C:\Users\devora.CHELEM\Sources\csm\test_cases\inbal\reading fragments'

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

#pdb without chains
#pdb with chains
#pdb with hetatm
#pdb with hetatm no flag
#pdb with use-sequence
def test_pdb_with_hetatm():
    file=os.path.join(inbal_folder, 'water-6.pdb')
    m=Molecule.from_file(file, format='pdb', read_fragments=True, use_chains=True)
    assert(len(m.chains)==6)

def test_pdb_with_hetatm_no_frag():
    file=os.path.join(inbal_folder, 'water-6.pdb')
    m=Molecule.from_file(file, format='pdb', use_chains=True)
    assert(len(m.chains)==6)

def test_pdb_with_hetatm_no_chain():
    file=os.path.join(inbal_folder, 'water-6.pdb')
    m=Molecule.from_file(file, format='pdb')
    assert(len(m.chains)==1)

def test_mol_with_multiple():
    file=os.path.join(inbal_folder, 'water-6.mol')
    m=Molecule.from_file(file, format='mol', read_fragments=True, use_chains=True)
    assert(len(m.chains)==6)

def test_mol_with_multiple_no_frag():
    file = os.path.join(inbal_folder, 'water-6.mol')
    m = Molecule.from_file(file, format='mol', use_chains=True)
    assert (len(m.chains) == 1)

def test_mol_with_multiple_no_chain():
    file = os.path.join(inbal_folder, 'water-6.mol')
    m = Molecule.from_file(file, format='mol')
    assert (len(m.chains) == 1)

#test for the following calculation types:
#approx with use-chains
#approx without use-chains
#trivial with use-chains
#trivial without use-chains
#exact with keep-structure
#normalizations with use-chains, 0-6