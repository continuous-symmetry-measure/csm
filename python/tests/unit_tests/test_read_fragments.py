import os
import pytest
from csm.molecule.molecule import Molecule, MoleculeReader
from csm.input_output.arguments_old import get_split_arguments
import logging
from conftest import test_folder



class baseClass:
    format="pdb"
    def no_flags(self):
        args = get_split_arguments(['c2', self.path, 'bla.txt'])
        m = MoleculeReader.from_file(**args)
        assert (len(m.chains) == 1)
        return m

    def read_fragments(self, caplog):
        args = get_split_arguments(['c2', self.path, 'bla.txt', '--read-fragments'])
        assert args['use_chains']==True
        m = MoleculeReader.from_file(**args)
        return m

    def use_chains(self, caplog):
        args = get_split_arguments(['c2', self.path, 'bla.txt', '--use-chains'])
        m = MoleculeReader.from_file(**args)
        return m

    def read_fragments_use_chains(self):
        args = get_split_arguments(['c2', self.path, 'bla.txt', '--read-fragments', '--use-chains'])
        m = MoleculeReader.from_file(**args)
        return m


class TestMolWithMultiple(baseClass):
    path = os.path.join(test_folder, 'water-6.mol')
    format = "mol"

    # scenario: mol file with 3 molecules separated by $$$,  user selected no flags
    # result: a single molecule matching the first molecule in the file
    def test_no_flags(self):
        m = super().no_flags()
        assert len(m) == 3

    # scenario: mol file with 3 molecules separated by $$$,  user selected --use-chains
    # result: a single molecule matching the first molecule in the file
    def test_use_chains(self, caplog):
        m = super().use_chains(caplog)
        assert len(m) == 3

    # scenario: mol file with 3 molecules separated by $$$,  user selected --read-fragments
    # result: a single molecule composed of 3 chains, each one matching the atoms described for their respective molecule in the file
    def test_read_fragments(self, caplog):
        m = super().read_fragments(caplog)
        assert len(m) == 18
        assert len(m.chains) == 6

    # scenario: mol file with 3 molecules separated by $$$,  user selected --read-fragments, --use-chains
    # result: a single molecule composed of 3 chains, each one matching the atoms described for their respective molecule in the file (the --use-chains is redundant, as it is assumed by --read-fragments)
    def test_read_fragments_use_chains(self):
        m = super().read_fragments_use_chains()
        assert len(m) == 18
        assert len(m.chains) == 6


class TestModelEndMdlWithIDS(baseClass):
    path=os.path.join(test_folder, 'model-endmdl-withIDs.pdb')
    def test_no_flags(self):
        m=super().no_flags()
        assert len(m)==6
        assert len(m.chains)==1

    #scenario: pdb with three models, each with 3 chains. user selected read-fragments
    #result: a single molecule with three chains matching the three models. No chains from the internal chains.
    def test_read_fragments(self, caplog):
        m=super().read_fragments(caplog)
        assert len(m)==18
        assert len(m.chains)==3

    # scenario: pdb with three models, each with 3 chains. user selected use-chains
    # result: the first molecule, with three chains matching the three internal chains.
    def test_use_chains(self, caplog):
        m=super().use_chains(caplog)
        assert len(m)==6
        assert len(m.chains)==3

    #scenario: pdb with three models, each with 3 chains. user selected read-fragments, use-chains
    #result: a single molecule with three chains matching the three models. No chains from the internal chains.
    def test_read_fragments_use_chains(self, caplog):
        m=super().read_fragments_use_chains()
        assert len(m)==18
        assert len(m.chains)==3


class TestHetAtmWithIDS(baseClass):
    path = os.path.join(test_folder, 'hetatm-with-IDs.pdb')

    def test_no_flags(self):
        m=super().no_flags()
        assert len(m)==18
        assert len(m.chains)==1

    def test_read_fragments(self, caplog):
        m=super().read_fragments(caplog)
        assert len(m)==18
        assert len(m.chains)==1

    def test_use_chains(self, caplog):
        m=super().use_chains(caplog)
        assert len(m)==18
        assert len(m.chains)==6

    def test_read_fragments_use_chains(self):
        m=super().read_fragments_use_chains()
        assert len(m)==18
        assert len(m.chains)==1


class TestMolWithSingle(baseClass):
    path=os.path.join(test_folder, "just-one-mol.mol")
    format="mol"

    def test_no_flags(self):
        m=super().no_flags()
        assert len(m)==3
        assert len(m.chains)==1

    def test_read_fragments(self, caplog):
        m=super().read_fragments(caplog)
        assert len(m)==3
        assert len(m.chains)==1

    def test_use_chains(self, caplog):
        m=super().use_chains(caplog)
        assert len(m)==3
        assert len(m.chains)==1

    #scenario: mol file with 1 molecule,  user selected --read-fragments, --use-chains
    #result: a single molecule matching the 1 molecule, no chains (and a message to user that file has no fragments)
    def test_read_fragments_use_chains(self):
        m=super().read_fragments_use_chains()
        assert len(m)==3
        assert len(m.chains)==1



class xTestModelEndMdlWithoutIDS(baseClass):
    path = os.path.join(test_folder, 'model-endmdl-withoutIDs.pdb')
    #This is supposedly a broken pdb and doesn't work, but for some reason it does, so I don't get what's going on
    def test_no_flags(self):
        m= super().no_flags()




class xTestPDBWithMixOFHETATMandATOM(baseClass):
    file=""
    #scenario: pdb model that contains chains and hetatms. user selected no flags
    #result: a single molecule with no chains, containing all the atoms marked ATOM and HETATM in the pdb file

    #scenario: pdb model that contains chains and hetatms. user selected --use-chains
    #result: a single molecule with chains matching those indicated in the pdb file, including the chains of the HETATMs

    #scenario: pdb model that contains chains and hetatms. user selected --use-chains, --read-fragments
    #result: a single molecule with no chains, and a mesage to user that the file has no fragments







class xTestTER(baseClass):
        path = os.path.join(test_folder, 'TER.pdb')

        def test_read_fragments(self, caplog):
            m = super().test_read_fragments(caplog)
            assert len(m) == 18
            assert len(m.chains) == 6



#scenario: file with hetatm and chainid
#with --use-chains, the hetatms will be assigned chains accordingly to chainID
# with --read-fragments, the chain-id will be ignored.

#file with hetatm without chainID and endmodel
#result: openbabel error

#file with hetatm and endmdl
#with --use-chains, the first model will be read, and hetatms will be assigned chains according to their chain id.
#with --read-fragments, all the molecules will be read as chains, and the ID will be ignored.

#file with mixture of alphabetic HETATM and ATOM
#The program expects Hetatm chain IDS to always be in column 26.
#If the Hetatm chain id is not in column 26, this will currently cause problems. (program probabaly won't crash, but rather will provide a wrong chain ID)

#file with mixture of non-alphabetic HETATM, and ATOM
#the program does not recognize this as a problem, and will assign the Hetatm a chainID from column 26 and continue onwards. It is the user's responsibility to fix the file before running csm.

#file with mixture of chainID and endmdl
#if the user selects --use-chains, only the first molecule will be read, and the chainID information will be used
#if the user selects --read-fragments, the molecules will be read as chains, and the internal chain information will be ignored.
