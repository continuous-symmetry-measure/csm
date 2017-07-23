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


#mol with fragments, no flags
def test_should_only_read_first_molecule_without_read_fragments():
    file = os.path.join(inbal_folder, 'water-6.mol')
    m = Molecule.from_file(file, format='mol')
    assert (len(m.chains) == 1)
    assert (len(m)==3)

#mol with fragments, --use-chains
#print: no chains
#How can i test if there is a printout?
def test_should_only_read_first_molecule_without_read_fragments_ignore_use_chains(caplog):
    file = os.path.join(inbal_folder, 'water-6.mol')
    m = Molecule.from_file(file, format='mol', use_chains=True)
    for record in caplog.records():
        hi=1
    assert (len(m.chains) == 1)
    assert (len(m)==3)

#mol with fragments, --read-fragments
#print: added use chains automatically
def test_should_read_molecules_as_chains_with_read_fragments(capsys):
    file=os.path.join(inbal_folder, 'water-6.mol')
    m=Molecule.from_file(file, format='mol', read_fragments=True)
    out, err = capsys.readouterr()
    assert(len(m.chains)==6)
    assert(len(m)==18)

#mol with fragments, --use-chains, -read-fragments
def test_should_read_molecules_as_chains_with_read_fragments():
    file=os.path.join(inbal_folder, 'water-6.mol')
    m=Molecule.from_file(file, format='mol', read_fragments=True, use_chains=True)
    assert(len(m.chains)==6)
    assert(len(m)==18)

#pdb with endmdl, no flags

#pdb with endmdl, --use-chains

#pdb with endmdl, --read-fragments

#pdb with endmdl, --use-chains, --read-fragments

#pdb with chains(ter), no flags

#pdb with chains(ter), --use-chains

#pdb with chains(ter), --read-fragments

#pdb with chains(ter), --use-chains, --read-fragments

#pdb with chains and endmdl, flags irrelevant:


#pdb with endmdl, --ignore-hetatm

#pdb with endmdl, --use-chains, --ignore-hetatm

#pdb with endmdl, --read-fragments, --ignore-hetatm

#pdb with endmdl, --use-chains, --read-fragments, --ignore-hetatm

#pdb with chains(ter), --ignore-hetatm

#pdb with chains(ter), --use-chains, --ignore-hetatm

#pdb with chains(ter), --read-fragments, --ignore-hetatm

#pdb with chains(ter), --use-chains, --read-fragments, --ignore-hetatm




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
    assert (len(m)==3)

def test_mol_with_multiple_no_chain():
    file = os.path.join(inbal_folder, 'water-6.mol')
    m = Molecule.from_file(file, format='mol')
    assert (len(m.chains) == 1)

def test_pdb_with_chains_hetatm_and_nonequal_equivalence_classes():
    file=r'D:\UserData\devora\Sources\csm\test_cases\old_test_cases\inbal\proteins\2xql.pdb'
    m=Molecule.from_file(file, format='pdb', use_chains=True)
#test for the following calculation types:
#approx with use-chains
#approx without use-chains
#trivial with use-chains
#trivial without use-chains
#exact with keep-structure
#normalizations with use-chains, 0-6

#add ignore_hetatms

#scenario: mol file with 3 molecules separated by $$$,  user selected no flags
#result: a single molecule matching the first molecule in the file

#scenario: mol file with 3 molecules separated by $$$,  user selected --use-chains
#result: a single molecule matching the first molecule in the file


#scenario: mol file with 3 molecules separated by $$$,  user selected --read-fragments
#result: a single molecule composed of 3 chains, each one matching the atoms described for their respective molecule in the file
#(currently this scenario raises an error and asks that the user also select --use-chains. since I understood from your email that that's not the desired behavior, I would change this so that inputting --read-fragments automatically selects --use-chains as well)
#(theoretical alternate result: a single molecule with the atoms from all three molecules and no chains. it sounds like you don't want this)

#scenario: mol file with 3 molecules separated by $$$,  user selected --read-fragments, --use-chains
#result: a single molecule composed of 3 chains, each one matching the atoms described for their respective molecule in the file

#scenario: mol file with 1 molecule,  user selected --read-fragments, --use-chains
#result: a single molecule matching the 1 molecule. no chains
#prints message

#scenario: pdb model that contains chains and hetatms. user selected no flags
#result: a single molecule with no chains, containing all the atoms marked ATOM and HETATM in the pdb file

#naming chains-- other than chains in pdb, all chains get a oridnal number

#scenario: pdb model that contains chains and hetatms. user selected --use-chains
#result: a single molecule with chains matching those indicated in the pdb file. should the hetatms also be marked as chains?

#scenario: pdb model that contains chains and hetatms. user selected --use-chains, --read-fragments
#result: a single molecule with chains matching those indicated in the pdb file, plus hetatms as chains??

#scenario: pdb with two models, no chains. user selected read-fragments
#result: a single molecule with two chains matching the two models


#scenario: pdb with two models, each with 2 chains. user selected read-fragments
#result: a single molecule with two chains matching the two models, and no chaisn from the chains


#scenario: pdb with two models, each with 2 chains. user selected read-fragments, use-chains
#result: a single molecule with 4 chains?