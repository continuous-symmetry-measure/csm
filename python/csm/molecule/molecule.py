"""
@author: Devora Witty
"""
import json
import logging
from collections import OrderedDict
from pathlib import Path

import copy
import numpy as np
import os
try:
    from openbabel.openbabel import OBAtomAtomIter, OBConversion, OBMol, OBMolAtomIter, obErrorLog, obError
except ImportError:
    from openbabel import OBAtomAtomIter, OBConversion, OBMol, OBMolAtomIter, obErrorLog, obError

from csm.input_output.formatters import csm_log as print
from csm.input_output.formatters import silent_print
from csm.molecule.atom import Atom, GetAtomicSymbol
from csm.molecule.normalizations import normalize_coords, de_normalize_coords, calculate_norm_factor

logger = logging.getLogger("csm")

ob_debug = False
if not ob_debug:
    obErrorLog.SetOutputLevel(obError)


def get_format(format, filename):
    if not format:
        format =  Path(filename).suffix.strip(".")
    if format.lower() == "csm":
        return "csm"
    conv = OBConversion()
    if not conv.SetInFormat(format):
        format=None
    if not format:
        raise ValueError("Error discovering format from filename " + filename)
    return format


class Chains(OrderedDict):
    '''
    two sets of keys, string keys and integer keys that match the index of the chain
    the inner dict is composed of integer keys, string keys must be translated to int before sending on to inner dict
    '''

    def __init__(self):
        self._indexes_to_strings = []
        self._strings_to_indexes = {}
        super().__init__()

    def __getitem__(self, key):
        if isinstance(key, int):
            return super().__getitem__(key)
        if isinstance(key, str):
            return super().__getitem__(self._strings_to_indexes[key])
        raise ValueError

    def __setitem__(self, key, value):
        '''
        :param key: expects a string and will force ints to string
        :param value:
        :return:
        '''
        if isinstance(key, int):
            key = str(key)
        self._indexes_to_strings.append(key)
        index = self._indexes_to_strings.index(key)
        self._strings_to_indexes[key] = index
        return super().__setitem__(index, value)

    def __contains__(self, item):
        if isinstance(item, int):
            return item < len(self._indexes_to_strings)
        if isinstance(item, str):
            try:
                return self._strings_to_indexes.__contains__(item)
            except KeyError:
                return False

    def to_array(self):
        zipped_arr = [(key, self.__getitem__(key)) for key in self._indexes_to_strings]
        return zipped_arr

    def from_array(self, arr_of_tuples):
        for key, val in arr_of_tuples:
            self.__setitem__(key, val)

    def index_to_name(self, index):
        return self._indexes_to_strings[index]

    def name_to_index(self, name):
        test= self._strings_to_indexes[name]
        return test


class MoleculeMetaData:
    '''
    This class is primarily used to store metadata needed to write results, although format+filecontent+babel_bond
    are sometimes used to recreate a molecule from scratch, see: redo_molecule

    format, filename, and babel_bond are set in _initialize_single_molecule
    file_content is set in read_obm or read_csm
    index is set in read_multiple_molecules, or, revoltingly, in do_commands after calling redo_molecule
    '''

    def __init__(self, file_content=[], format=None, filepath="", babel_bond=False, index=0, initial_title="",
                 initial_comments="", use_filename=True, out_format=None, selected_mols=[]):
        self.file_content = file_content
        self.format = format
        self._out_format = out_format
        self.filepath = filepath
        self.babel_bond = babel_bond
        self.index = index
        self.initial_title = initial_title
        self.initial_comments = initial_comments
        self.use_filename = use_filename
        self.select_mols = selected_mols

    @property
    def filename(self):
        return os.path.basename(self.filepath)

    @property
    def out_format(self):
        if self._out_format:
            return self._out_format
        return self.format

    @staticmethod
    def from_dict(self, dict):
        m = MoleculeMetaData(**dict)
        return m

    def to_dict(self):
        return vars(self)

    def appellation(self, no_file_format=False, no_leading_zeros=False, write_original_mols_index=False):
        '''
        the name for the molecule when printing to screen or creating tables.
        if the molecule was read from a folder, it's the filename
        otherwise it's the internal index
        :param no_file_format: remove file ending (eg .xyz)
        :param no_leading_zeros: print the index without leading zeroes
        :param write_original_mols_index: boolean value for use the original indexes for the selected molecules.
        :return: 
        '''
        if self.use_filename:
            if no_file_format:
                return Path(self.filename).stem
            return self.filename

        if write_original_mols_index:
            mol_index = self.select_mols[self.index] + 1   # start from 1 instead of 0
        else:
            mol_index = self.index + 1  # start from 1 instead of 0
        if no_leading_zeros:
            return str(mol_index)

        mol_str = "%04d" % mol_index
        return mol_str


class Molecule:
    """
    Represents a molecule for CSM calculation.
    """

    def __init__(self, atoms=[], to_copy=False):
        """  
        :param atoms: a list of Atoms, default empty
        :param to_copy: boolean, default False, when True none of the molecules fields will be filled
        """
        if not to_copy:
            if len(atoms)<1: raise ValueError("Cannot create molecule with no atoms")
            self._atoms = atoms
            self._norm_factor = 1.0  # is recalculated by molecule._complete_initialization

            self._calculate_center_of_mass()
            self._create_bondset()
            self.create_Q()

            # not included:
            # equivalence class initialization, chain initialization
            self._equivalence_classes = []
            self._chains = Chains()
            self._groups_with_internal_chains = []
            self._chains_with_internal_groups = {}
            self._chain_equivalences = []

            # getting rid of obmol:
            # self._obmol = obmol
            self._deleted_atom_indices = []
            self.has_been_normalized = None

            self.metadata = MoleculeMetaData()

    def __str__(self):
        return self.metadata.appellation()

    def copy(self):
        # deepcopy is used only for atoms,
        # because atoms are the only property changed between runs of Directions that necessitated copying
        # (due to call to denormalize)

        m = Molecule(to_copy=True)
        m._atoms = copy.deepcopy(self.atoms)
        m._norm_factor = self.norm_factor
        m.has_been_normalized = self.has_been_normalized

        m._deleted_atom_indices = self._deleted_atom_indices
        m.metadata = self.metadata

        m._bondset = self.bondset
        m._Q = self.Q
        m._center_of_mass = self._center_of_mass

        m._equivalence_classes = self.equivalence_classes
        m._chains = self._chains
        m._groups_with_internal_chains = self._groups_with_internal_chains
        m._chains_with_internal_groups = self._chains_with_internal_groups
        m._chain_equivalences = self._chain_equivalences

        return m

    def to_dict(self):
        return {
            # critical to include
            "atoms": [atom.to_dict() for atom in self.atoms],

            # trivial to include
            "norm factor": self.norm_factor,
            "center of mass": self.center_of_mass,
            "normalized": self.has_been_normalized,

            # expensive to recalculate
            "equivalence classes": self.equivalence_classes,  # the most expensive part of loading large molecule
            "groups_with_internal_chains": self.groups_with_internal_chains,
            "chains_with_internal_groups": self.chains_with_internal_groups,
            "chain_equivalences": self.chain_equivalences,

            # Classes:
            # obmol: needed for printing:
            "deleted indices": self._deleted_atom_indices,
            "metadata": self.metadata.to_dict(),

            # chains
            "chains": self.chains.to_array()

            # unsure whether worth the bother of including
            # "bondset":self.bondset,
            # "Q":self.Q, #almost definitely not worth the bother
        }

    @staticmethod
    def from_json(raw_json):
        json_dict = json.loads(raw_json)
        return Molecule.from_dict(json_dict)

    @staticmethod
    def from_dict(in_dict):
        atoms = [Atom.from_dict(a) for a in in_dict["atoms"]]
        m = Molecule(atoms)

        m.norm_factor = in_dict["norm factor"]
        m.has_been_normalized = in_dict["normalized"]

        c = Chains()
        c.from_array(in_dict["chains"])
        m._chains = c

        m._deleted_atom_indices = in_dict["deleted indices"]
        m.metadata = MoleculeMetaData.from_dict(in_dict["metadata"])

        m._center_of_mass = in_dict["center of mass"]
        m._equivalence_classes = in_dict["equivalence classes"]
        unfixed_groups = in_dict["groups_with_internal_chains"]
        m._groups_with_internal_chains = []
        for group in unfixed_groups:
            fixed = {int(key): value for key, value in group.items()}
            m._groups_with_internal_chains.append(fixed)
        unfixed_chains = in_dict["chains_with_internal_groups"]
        fixed = {int(key): value for key, value in unfixed_chains.items()}
        m._chains_with_internal_groups = fixed
        m._chain_equivalences = in_dict["chain_equivalences"]

        m.create_Q()
        m._create_bondset()
        return m

    @property
    def center_of_mass(self):
        return self._center_of_mass

    @property
    def Q(self):
        '''
        :return: a numpy array of the molecule's atoms' coordinates, each represented as a numpy array
        '''
        return self._Q

    @property
    def atoms(self):
        '''
        :return: a list of atom objects
        '''
        return self._atoms

    @property
    def bondset(self):
        '''
        :return: a set of tuples, each representing a pair of atoms with a bond between them
        '''
        return self._bondset

    @property
    def equivalence_classes(self):
        '''
        :return: a list of the molecule's equivalence classes (themselves represented as lists), sorted in order of length
        '''
        self._equivalence_classes.sort(key=len)
        return self._equivalence_classes

    @property
    def norm_factor(self):
        # Normalization factor. Defaults to 1.0 if normalize wasn't called
        return self._norm_factor

    @norm_factor.setter
    def norm_factor(self, nf):
        self._norm_factor = nf

    @property
    def chains(self):
        '''
        :return: a dictionary of the molecule's chains, with keys being the name of the chains
        '''
        return self._chains

    @property
    def groups_with_internal_chains(self):
        return self._groups_with_internal_chains

    @property
    def chains_with_internal_groups(self):
        return self._chains_with_internal_groups

    @property
    def chain_equivalences(self):
        return self._chain_equivalences

    def __len__(self):
        return len(self._atoms)

    def _create_bondset(self):
        self._bondset = set()
        for i in range(len(self._atoms)):
            for match in self._atoms[i].adjacent:
                self._bondset.add((i, match))

    def _initialize_chains(self, use_chains):
        """
        1. Create chains, and overwrite atom.chain with index of chain rather than string
        2. within each equivalence class, labels by chain
            self._groups_with_internal_chains=[array of equivalence classes:
                            {dictionary of chains:
                                [array of indexes belonging to that chain in that equivalence class]
                            }
                        ]
            self._chains_with_internal_groups= {dicitonary of chains:
        3. Calculates chain equivalencies. (e.g, chains A and B are equivalent, chains C and D are equivalent)
        """

        # create chains, and overwrite atom.chain with index of chain rather than string
        chains = Chains()
        for i, atom in enumerate(self.atoms):
            chain = atom.chain
            if not use_chains:
                chain = 'Simulated chain'
            if chain not in chains:
                chains[chain] = []
            chains[chain].append(i)
            atom.chain = chains._strings_to_indexes[chain]
        self._chains = chains

        # within each equivalence class, labels by chain,
        # within each chain, group by equivalence class
        groups_with_internal_chains = []
        num_equiv = len(self.equivalence_classes)
        chains_with_internal_groups = {chain: [None] * num_equiv for chain in self.chains}

        for group_index, group in enumerate(self.equivalence_classes):
            chaingroup = {}
            for atom_index in group:
                chain_index = self.atoms[atom_index].chain

                # add index to array in chaingroup dict
                try:
                    chaingroup[chain_index].append(atom_index)
                except KeyError:
                    chaingroup[chain_index] = [atom_index]
                # add index to classed_chains
                try:
                    chains_with_internal_groups[chain_index][group_index].append(atom_index)
                except AttributeError:
                    chains_with_internal_groups[chain_index][group_index] = [atom_index]
            # add the completed chaingroup dict to group_chains arrat
            groups_with_internal_chains.append(chaingroup)

        self._groups_with_internal_chains = groups_with_internal_chains
        self._chains_with_internal_groups = chains_with_internal_groups

        # Calculates chain equivalencies. (e.g, chains A and B are equivalent, chains C and D are equivalent)
        # we define chains as equivalent if the number of atoms they have in each equivalence class are the same
        chain_equivalences = []
        marked = [False] * len(self.chains)

        for chain_index, chainkey in enumerate(self.chains):
            if marked[chain_index]:
                continue
            equiv = []
            for chain2_index, chainkey2 in enumerate(self.chains):
                # start with a simple length check to spare checking equivalence classes if chains arent same length to begin with
                # if len(self.chains[chainkey])!=len(self.chains[chainkey2]):
                #    continue
                same_lengths = True
                for group in groups_with_internal_chains:
                    if not same_lengths:
                        break
                    try:
                        length = len(group[chain_index])
                    except KeyError:  # that chain is notn in this equivalence group
                        continue
                    try:
                        if length != len(group[chain2_index]):
                            same_lengths = False
                    except KeyError:
                        same_lengths = False
                if same_lengths:
                    equiv.append(chain2_index)
                    marked[chain2_index] = True
            chain_equivalences.append(equiv)
        self._chain_equivalences = chain_equivalences

    def _calculate_equivalency(self):
        """
        Preprocess a molecule based on the arguments passed to CSM
        """

        def is_similar(atoms_group_num, a, b):
            found = True
            mark = set()

            valency_a = len(self._atoms[a].adjacent)
            valency_b = len(self._atoms[b].adjacent)

            # for each of i's neighbours
            for i in range(valency_a):
                found = False

                for j in range(valency_b):
                    if j in mark:
                        continue

                    if atoms_group_num[self._atoms[a].adjacent[i]] == atoms_group_num[self._atoms[b].adjacent[j]]:
                        # the i-th neighbour of 'a' belongs to the same group as the j-th neighbour of 'b'
                        found = True
                        mark.add(j)
                        break

                if not found:
                    break

            return found

        group_num = 0
        groups = []
        atoms_size = len(self._atoms)
        marked = set()

        atoms_group_num = {}

        # break into initial groups by symbol and valency
        for i in range(atoms_size):
            if i in marked:
                continue

            groups.append([])

            for j in range(atoms_size):
                if j in marked \
                        or len(self._atoms[i].adjacent) != len(self._atoms[j].adjacent) \
                        or self._atoms[i].symbol != self._atoms[j].symbol:
                    continue

                groups[group_num].append(j)
                atoms_group_num[j] = group_num
                marked.add(j)

            group_num += 1

        # logger.debug("initial number of groups:" + str(group_num))
        # iteratively refine the breakdown into groups
        # break into subgroups at an infinite depth - as long as there's something to break, it is broken

        divided_group = True
        num_iters = 0

        while divided_group:
            num_iters += 1
            divided_group = False

            new_groups = []

            for i, group in enumerate(groups):
                first_elem = group[0]
                group_size = len(group)

                sub_group = []

                # for each item in the group (except the first) check if it can be split or not
                for j in range(1, group_size):
                    if not is_similar(atoms_group_num, group[j], first_elem):
                        # add elem to new subGroup
                        sub_group.append(group[j])
                        group[j] = -1

                if len(sub_group) > 0:
                    divided_group = True
                    new_groups.append(sub_group)
                    for el in sub_group:
                        atoms_group_num[el] = group_num
                    group_num += 1
                    # remove elements of sub_group from group
                    groups[i] = [el for el in group if el != -1]

            groups.extend(new_groups)

        # logger.debug("Broken into %d groups with %d iterations." % (group_num, num_iters))

        self._equivalence_classes = groups
        for group in groups:
            for atom_index in group:
                for equiv_index in group:
                    self._atoms[atom_index].add_equivalence(equiv_index)

    def strip_atoms(self, remove_hy=False, select_atoms=[], ignore_atoms=[], use_backbone = False):
        """
            Creates a new Molecule from m by removing atoms who's symbol is in the remove list
            :param csm_args:
            :param removeList: atomic symbols to remove
        """

        indices_to_remove=[]
        if select_atoms:
            indices_to_remove=[i for i in range(len(self._atoms)) if i not in select_atoms]
        elif ignore_atoms:
            indices_to_remove = ignore_atoms
        if use_backbone:
            backbone_atoms = ['N', 'CA', 'C', 'O']
            indices_to_remove.extend([i for i in range(len(self._atoms)) if self._atoms[i].atom_name not in backbone_atoms])
        indices_to_remove=set(indices_to_remove)

        #check for bad input 1: index provided that doesnt exist:
        set_of_all_atoms = set(range(len(self._atoms)))
        if len(set_of_all_atoms.union(indices_to_remove))!=len(set_of_all_atoms):
            raise ValueError("An atom index you have input to --select-atoms or --ignore-atoms does not exist in the molecule")

        #add remove hy:
        if remove_hy:
            hy_atoms_indices=set([i for i in range(len(self._atoms)) if self._atoms[i].symbol =="H"])
            #add the hydrogens
            indices_to_remove=indices_to_remove.union(hy_atoms_indices)

        indices_to_remove=sorted(indices_to_remove)

        self._deleted_atom_indices=indices_to_remove #needed when writing output of molecule to file via openbabel
        self.fixed_indexes = fixed_indexes= [i for i in range(len(self))] #initializing here before the return because we use it when reading permutation
        if not indices_to_remove: #relevant if remove-hy was selected but no hydrogen in molecule, or select-atoms selcted every atom
            return #(may as well save time)

        #we now store a mapping of each atom to its new index once all relevant atoms have been removed, in order to update connectivity
        num_removed_atoms=0
        for i in range(len(self._atoms)):
            if i in indices_to_remove:
                num_removed_atoms+=1
                fixed_indexes[i] = None
            else:
                # however many atoms have been removed up to this index is the amount its index needs adjusting by
                fixed_indexes[i] -= num_removed_atoms


        # adjust the connectivity indices before we do any popping whatsoever
        for i, atom in enumerate(self._atoms):
            adjacent_new = []
            for adjacent in atom.adjacent:
                if fixed_indexes[adjacent] is not None:
                    adjacent_new.append(fixed_indexes[adjacent])
            atom.adjacent = adjacent_new

        for to_remove in reversed(indices_to_remove):  # reversed order because popping changes indexes after
            self._atoms.pop(to_remove)

        self.fixed_indexes = fixed_indexes #overwrite default value
        self._create_bondset()

        # logger.debug(len(removed_atoms), "atoms removed")

    def _calculate_center_of_mass(self):
        coords = [atom.pos for atom in self._atoms]
        masses = [atom.mass for atom in self._atoms]
        x_avg = y_avg = z_avg = 0.0

        size = len(masses)

        mass_sum = 0
        for i in range(size):
            x_avg += coords[i][0] * masses[i]
            y_avg += coords[i][1] * masses[i]
            z_avg += coords[i][2] * masses[i]
            mass_sum += masses[i]
        x_avg /= mass_sum
        y_avg /= mass_sum
        z_avg /= mass_sum

        self._center_of_mass = [x_avg, y_avg, z_avg]

    def normalize(self):
        """
        Normalize the molecule
        """
        coords = [atom.pos for atom in self._atoms]
        self._calculate_center_of_mass()
        self._norm_factor = calculate_norm_factor(coords, self.center_of_mass)
        norm_coords = normalize_coords(coords, self.center_of_mass, self.norm_factor)
        size = len(self._atoms)
        for i in range(size):
            self._atoms[i].pos = norm_coords[i]
        self.create_Q()
        self.has_been_normalized = True

    def de_normalize(self):
        coords = [atom.pos for atom in self._atoms]
        denorm_coords = de_normalize_coords(coords, self.norm_factor)

        size = len(self._atoms)
        for i in range(size):
            self._atoms[i].pos = denorm_coords[i]
        self.create_Q()
        self.has_been_normalized = False

    def create_Q(self):
        self._Q = np.array([np.array(atom.pos) for atom in self.atoms])

    def print_equivalence_class_summary(self, display_chains):
        """
        Displays information about equivalence classes and chains
        """
        lengths = {}
        for group in self._equivalence_classes:
            try:
                lengths[len(group)] += 1
            except KeyError:
                lengths[len(group)] = 1
        for key in lengths:
            silent_print(
                "%d group%s of length %d" % (lengths[key], 's' if lengths[key] and lengths[key] > 1 else '', key))

        if display_chains:
            if len(self.chains) > 1:
                for chain in self.chains:
                    silent_print(
                        "Chain %s of length %d" % (self.chains._indexes_to_strings[chain], len(self.chains[chain])))
                silent_print(
                    "%d equivalence class%s of chains" % (len(self.chain_equivalences), 'es' if lengths[key] else ''))
                for chaingroup in self.chain_equivalences:
                    chainstring = "Group of length " + str(len(chaingroup)) + ":"
                    for index in chaingroup:
                        chainstring += " "
                        chainstring += str(self.chains._indexes_to_strings[index])
                    silent_print(str(chainstring))
            # else:
            #    silent_print("Molecule has no chains")

    def _complete_initialization(self, use_chains, remove_hy, select_atoms=[], ignore_atoms=[], use_backbone=False):
        """
        Finish creating the molecule after reading the raw data
        """
        self.strip_atoms(remove_hy, select_atoms, ignore_atoms, use_backbone)
        self._calculate_equivalency()
        self._initialize_chains(use_chains)
        self.normalize()

    @staticmethod
    def xyz_string(atoms, positions=None, header=""):
        '''
        :param atoms: an array of atoms
        :param positions: optional, an array of positions to overwrite positions in the array of atoms
        :return:
        '''
        coords = str(len(atoms))
        coords += "\n" + header + "\n"
        for i, atom in enumerate(atoms):
            coords += str(atom.symbol)
            for index, coor in enumerate(atom.pos):
                coords += "\t"
                if positions is not None:
                    coords += str(positions[i][index])
                else:
                    coords += str(coor)
            coords += "\n"
        return coords


class MoleculeFactory:
    @staticmethod
    def dummy_molecule_from_size(size, groups):
        '''
        The dummy molecule is a fake molecule used in the approx algorithm to generate chain permutations
        :param size: number of "atoms" in the fake molecule
        :param groups: the equivalence classes of those "atoms"
        :return:
        '''
        atoms = []

        for i in range(size):
            atom = Atom("C", (0, 0, 0), i, False)
            atoms.append(atom)
        mol = Molecule(atoms)
        mol._equivalence_classes = groups
        return mol

    @staticmethod
    def dummy_molecule_from_coords(coords, groups=None):
        atoms = []

        if groups is None:
            groups = [[i for i in range(len(coords))]]

        for fake_atom, group in enumerate(groups):
            for atom_index in group:
                coord = coords[atom_index]
                atom = Atom(str(fake_atom), coord, atom_index, False)
                atoms.append(atom)

        mol = Molecule(atoms)
        mol._calculate_equivalency()
        mol.normalize()
        return mol


class PDBLine:
    '''
    Currently only supports ATOM, HETATM, and CONECT records-- if any others are needed in the future they can be added then
    '''

    def __init__(self, pdb_line):
        self.pdb_line = pdb_line
        self._record_name = pdb_line[0:6]
        self.record_name = self._record_name.strip()

        self._atom_serial_number = pdb_line[6:11]  # used in read pdb connectivity
        if self.record_name in ["ATOM", "HETATM"]:
            self._init_atom_record(pdb_line)

        if self.record_name == "CONECT":
            self._init_conect_record(pdb_line)

    def _init_atom_record(self, pdb_line):
        self.full_atom_name = pdb_line[12:16]
        self.alternate_location_indicator = pdb_line[16]
        self.residue_name = pdb_line[17:20]

        # handle chains
        # ATOMxxxxxx1xxNxxxLYSxA
        # HETATMxxxx3xxHxxxHOHxxxxx1
        if self.record_name == 'ATOM':
            chain_designation = self.pdb_line[21]
        if self.record_name == 'HETATM':
            if self.pdb_line[21] != " ":
                chain_designation = self.pdb_line[21]
            else:
                chain_designation = self.pdb_line[25]
        self.chain_id = chain_designation

        self.sequence_number = pdb_line[22:26],  # used in use sequence
        self.iCode = pdb_line[26]
        self.x = pdb_line[30:38]
        self.y = pdb_line[38:46]
        self.z = pdb_line[46:54]
        self.occ = pdb_line[54:60]
        self.temp_factor = pdb_line[60:66]
        self.element = pdb_line[76:78]
        self.charge = pdb_line[78:80]

        self._parse_atom_name(pdb_line)

    def _parse_atom_name(self, pdb_line):
        # https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
        # Columns 13-16 - the atom name. Columns 13-14 contain the atom symbol, right aligned (except when the atom is for H).
        # Column 15 contains the remoteness indicator that represents the distance of the atom from the backbone of the protein. Greek letter remoteness codes are transliterated as follows:
        # alpha=A, beta=B, gamma = G, delta = D, epsilon=E, zeta= Z, eta=H].
        # Column 16 – is the branch designator
        # Column 17 - Alternate location.
        # *Atom names start with element symbols right-justified in columns 13-14 as permitted by the length of the name.
        # For example, the symbol FE for iron appears in columns 13-14, whereas the symbol C for carbon appears in column 14 (see Misaligned Atom Names).
        # If an atom name has four characters, however, it must start in column 13 even if the element symbol is a single character (for example, see Hydrogen Atoms).
        # Hydrogen atom records follow the records of all other atoms of a particular residue.
        # A hydrogen atom name starts with H. The next part of the name is based on the name of the connected nonhydrogen atom.
        # For example, in amino acid residues, H is followed by the remoteness indicator (if any) of the connected atom, followed by the branch indicator (if any) of the connected atom;
        # if more than one hydrogen is connected to the same atom, an additional digit is appended so that each hydrogen atom will have a unique name.
        # Hydrogen atoms in standard nucleotides and amino acids (other than the rarely seen HXT) are named according to the IUPAC recommendations (Pure Appl Chem 70:117 (1998) [abstract] [PDF]).
        # Names of hydrogen atoms in HETATM residues are determined in a similar fashion.
        # If the name of a hydrogen has four characters, it is left-justified starting in column 13; if it has fewer than four characters, it is left-justified starting in column 14.

        atom_name = pdb_line[12:16]
        atom_name = atom_name.strip()
        if atom_name[0] == "H":
            self.atom_symbol = "H"
            self.remoteness = ""
            self.branch_designation = ""
            self.duplicate_index = ""
            self.alternate_location = pdb_line[16]
            try:
                self.remoteness = self.atom_symbol[1]
                self.branch_designation = self.atom_symbol[2]
                self.duplicate_index = self.atom_symbol[3]
            except:  # TODO: check that what information is available gets recorded
                pass  # because it may be attempting to access information that doesn't exist so it's ok...
        else:
            self.atom_symbol = pdb_line[12:14].strip()
            self.remoteness = pdb_line[14]
            self.branch_designation = pdb_line[15]
            self.alternate_location = pdb_line[16]
        self.atom_name = pdb_line[12:15].strip()

        if not self.atom_symbol:
            raise ValueError(
                "PDB missing atom symbol information in atom name (columns 13-16). Please check that your pdb is formatted correctly")

    def _init_conect_record(self, pdb_line):
        adjacent_atoms = []
        stripped_pdb_line = pdb_line.strip()
        for i in range(11, len(stripped_pdb_line), 5):
            adjacent_atom_index = int(stripped_pdb_line[i:i + 5])
            adjacent_atoms.append(adjacent_atom_index)
        self.adjacent_atoms = adjacent_atoms

    @property
    def atom_serial_number(self):
        try:
            return int(self._atom_serial_number)
        except (ValueError, IndexError):
            return None

    @staticmethod
    def _pdb_line_to_dict(pdb_line):
        return PDBLine(pdb_line)


class MoleculeReader:
    """
    A static class that creates instances of Molecule from files or strings.
    """

    @staticmethod
    def _remove_multi_bonds(bonds):
        """
        Receives a list, returns a sorted list with no duplicates.
        :param bonds: a list
        :return: a sorted list with no duplicates
        """
        l = list(set(bonds))
        l.sort()
        return l

    @staticmethod
    def from_string(string, format, initialize=True,  # Never used, Chana
                    use_chains=False, babel_bond=False,
                    remove_hy=False, ignore_sym=False, use_mass=False):
        """
        create a Molecule from a string
        :param string: the string to create the molecule from
        :param format: the format of the string (any BabelBond supported format, eg "mol", "xyz")
        :param initialize: boolean, default True, when True equivalence classes and chains are calculated for the molecule
        :param use_chains: boolean, default False, when True chains are read from the string 
        :param babel_bond: boolean, default False, when True OpenBabel will attempt to guess connectivity information
        :param remove_hy: boolean, default False, when True hydrogen atoms will be removed
        :param ignore_sym: boolean, default False, when True all atom symbols will be "X" instead of elements 
                            (this affects equivalence class calculation)
        :param use_mass: boolean, default False, when True atomic mass will be used, when False all atoms have a mass of 1
        :return: an instance of class Molecule
        """
        # note: useMass is used when creating molecule, even though it is actually about creating the normalization
        # step one: get the molecule object
        obm = MoleculeReader._obm_from_strings([string], format, babel_bond)
        mol = MoleculeReader.mol_from_obm(obm, format, ignore_sym=ignore_sym, use_mass=use_mass)
        if initialize:
            mol._complete_initialization(use_chains, remove_hy)
        return mol

    @staticmethod
    def from_file(in_file_name, in_format=None, initialize=True,
                  use_chains=False, babel_bond=False,
                  remove_hy=False, ignore_sym=False, use_mass=False,
                  read_fragments=False, use_sequence=False,
                  keep_structure=False, select_atoms=[], conn_file=None,
                  out_format=None, ignore_atoms = [], use_backbone=False,
                  *args, **kwargs):
        """
        :param in_file_name: the name of the file to read the molecule from
        :param in_format: the format of the string (any BabelBond supported format, eg "mol", "xyz")
        :param initialize: boolean, default True, when True equivalence classes and chains are calculated for the molecule
        :param use_chains: boolean, default False, when True chains are read from the string 
        :param babel_bond: boolean, default False, when True OpenBabel will attempt to guess connectivity information
        :param remove_hy: boolean, default False, when True hydrogen atoms will be removed
        :param ignore_sym: boolean, default False, when True all atom symbols will be "X" instead of elements 
                            (this affects equivalence class calculation)
        :param use_mass: boolean, default False, when True atomic mass will be used, when False all atoms have a mass of 1
        :param read_fragments: boolean, default False, when True, multiple molecules in one file will be treated 
                                as "chains" of one composite molecule
        :param use_sequence: boolean, default False, for pdb files only-- uses pdb sequence information to establish 
                                equivalence (overwrites default equivalence)
        :param keep_structure: boolean, default False. If keep-structure is specified and molecule has no bond information, 
                                the peogram will raise an error and exit
        :return: an instance of class Molecule 
        """
        # suppress warnings from openbabel

        format = get_format(in_format, in_file_name)
        if format == "csm":
            mol = MoleculeReader._read_csm_file(in_file_name, ignore_sym, use_mass)
        else:
            obm = MoleculeReader._obm_from_file(in_file_name, format, babel_bond)
            mol = MoleculeReader.mol_from_obm(obm, format, ignore_sym=ignore_sym, use_mass=use_mass,
                                              read_fragments=read_fragments)
        mol = MoleculeReader._process_single_molecule(mol, in_file_name, format, initialize,
                                                      use_chains, babel_bond,
                                                      remove_hy, ignore_sym, use_mass,
                                                      read_fragments, use_sequence,
                                                      keep_structure, select_atoms, conn_file,
                                                      out_format, ignore_atoms, use_backbone)
        if not mol.bondset:
            if keep_structure:
                raise ValueError(
                    "User input --keep-structure but input molecule has no bonds. Did you forget --babel-bond?")
            else:
                print("Warning: Input molecule has no bond information")
        return mol

    @staticmethod
    def _process_single_molecule(mol, in_file_name, format, initialize=True,
                                 use_chains=False, babel_bond=False,
                                 remove_hy=False, ignore_sym=False, use_mass=False,
                                 read_fragments=False, use_sequence=False,
                                 keep_structure=False, select_atoms=[], conn_file=None,
                                 out_format=None, ignore_atoms=[], use_backbone=False, **kwargs):

        mol.metadata.format = format
        if out_format:
            mol.metadata._out_format = out_format
        mol.metadata.babel_bond = babel_bond
        mol.metadata.filepath = in_file_name

        if use_backbone and format.lower() != 'pdb':
            raise ValueError("--use-backbone only works with PDB files")

        if use_sequence:
            if format.lower() != 'pdb':
                raise ValueError("--use-sequence only works with PDB files")

            mol = MoleculeReader._create_pdb_with_sequence(mol, in_file_name, use_chains=use_chains,
                                                           babel_bond=babel_bond,
                                                           read_fragments=read_fragments, remove_hy=remove_hy,
                                                           ignore_sym=ignore_sym, use_mass=use_mass,
                                                           select_atoms=select_atoms, ignore_atoms=ignore_atoms,
                                                           use_backbone=use_backbone)
            # we initialize mol from within pdb_with_sequence because otherwise equivalnce classes would be overwritten
            return mol

        if format == "pdb":
            mol = MoleculeReader._read_pdb_connectivity_and_chains(in_file_name, mol, read_fragments, babel_bond, use_backbone)
        if conn_file and format == "xyz":
            MoleculeReader.read_xyz_connectivity(mol, conn_file)
        if initialize:
            mol._complete_initialization(use_chains, remove_hy, select_atoms, ignore_atoms, use_backbone)
            if len(mol.chains) < 2:
                if read_fragments:
                    print("Warning: Although you input --read-fragments, no fragments were found in file. "
                          "Fragments are marked by $$$ in mol files or by model/endmdl in pdb files")
                elif use_chains:
                    print("Warning: You specified --use-chains but molecule only has one chain")
        return mol

    @staticmethod
    def multiple_from_file(in_file_name, in_format=None, initialize=True,
                           use_chains=False, babel_bond=False,
                           remove_hy=False, ignore_sym=False, use_mass=False,
                           read_fragments=False, use_sequence=False,
                           keep_structure=False, select_atoms=[], conn_file=None,
                           out_format=None, ignore_atoms=[], use_backbone=False,
                           *args, **kwargs):
        #although the name of this function is "multiple from file", it is used both for files with multiple molecules
        #and for files with only a single molecule. it is used anytime the --input is a file, not a folder
        #it is extremely similar to .from_file. the difference is mostly in metadata, like molecule indices
        mols = []
        format = get_format(in_format, in_file_name)

        if format == "csm":
            mol = MoleculeReader._read_csm_file(in_file_name, ignore_sym, use_mass)
            mols.append(mol)

        else:
            obms = MoleculeReader._obm_from_file(in_file_name, format, babel_bond)
            if read_fragments:
                mol = MoleculeReader.mol_from_obm(obms, format, babel_bond=babel_bond, ignore_sym=ignore_sym,
                                                  use_mass=use_mass, read_fragments=read_fragments)
                if mol:  # stupid workaround for openbabel 3
                    mols.append(mol)

            else:
                obms, selected_mols=select_mols(obms, kwargs) #save a little bit of time- only continue processing the molecules that were selected
                for i, obm in enumerate(obms):
                    mol = MoleculeReader.mol_from_obm([obm], format, babel_bond=babel_bond, ignore_sym=ignore_sym,
                                                      use_mass=use_mass)
                    if mol: #stupid workaround for openbabel 3
                        mols.append(mol)

        processed_mols = []
        use_filename = True
        if len(mols) > 1 or kwargs.get("legacy_output"):
            use_filename = False
        for index, mol in enumerate(mols):
            p_mol = MoleculeReader._process_single_molecule(mol, in_file_name, format, initialize,
                                                            use_chains, babel_bond,
                                                            remove_hy, ignore_sym, use_mass,
                                                            read_fragments, use_sequence,
                                                            keep_structure, select_atoms, conn_file,
                                                            out_format, ignore_atoms,use_backbone)
            p_mol.metadata.index = index
            p_mol.metadata.use_filename = use_filename
            if not format == "csm" and not read_fragments:
                p_mol.metadata.select_mols=selected_mols
            processed_mols.append(p_mol)
        if mols and not p_mol.bondset:  # if mols -> checks if the mols list doesn't empty
            # this only checks for the final one,
            # on the assumption that all molecules in the file have the same bond status and to avoid printing a million times
            if keep_structure:
                raise ValueError(
                    "User input --keep-structure but input molecules have no bonds. Did you forget --babel-bond?")
            else:
                print("Warning: Input molecules have no bond information")

        return processed_mols

    @staticmethod
    def _obm_from_strings(strings, format, babel_bond=False):
        '''
        :param strings: an array of molecule strings
        :param format:
        :param babel_bond:
        :return:
        '''
        conv = OBConversion()
        obmol = OBMol()
        if not conv.SetInFormat(format):
            raise ValueError("Error setting openbabel format to" + format)
        if not babel_bond:
            conv.SetOptions("b", conv.INOPTIONS)
        obmols = []
        for string in strings:
            conv.ReadString(obmol, string)
            obmols.append(obmol)
            obmol = OBMol()
        return obmols

    @staticmethod
    def _obm_from_file(filename, format=None, babel_bond=False):
        """
        :param filename: name of file to open
        :param format: molecule format of file (eg xyz, pdb)
        :param babelBond:
        :return:
        """
        obmol = OBMol()
        conv = OBConversion()
        if not format:
            format = get_format(format, filename)
        if not conv.SetInFormat(format):
            raise ValueError("Error setting openbabel format to" + format)
        if format == "txt":
            raise ValueError("CSM does not support .txt format openbabel conversions")
        if not babel_bond:
            conv.SetOptions("b", conv.INOPTIONS)
        notatend = conv.ReadFile(obmol, filename)
        if not notatend:
            raise ValueError("Error reading file " + filename + " using OpenBabel, with format:", format)
        obmols = []
        while notatend:
            obmols.append(obmol)
            obmol = OBMol()
            notatend = conv.Read(obmol)
        return obmols

    @staticmethod
    def mol_from_obm(obmols, format, babel_bond=False, ignore_sym=False, use_mass=False, read_fragments=False,
                     **kwargs):
        """
        :param obmol: OBmol molecule
        :param args_dict: dictionary of processed command line arguments
        :return: A list of Atoms and a list of chains
        """
        if not read_fragments:
            obmols = obmols[:1]
        atoms = []
        mol_contents = []
        for obmol_id, obmol in enumerate(obmols):
            title_contents = obmol.GetTitle()
            mol_contents.append(mol_string_from_obm(obmol, format))
            for i, obatom in enumerate(OBMolAtomIter(obmol)):
                position = (obatom.GetX(), obatom.GetY(), obatom.GetZ())
                if ignore_sym:
                    symbol = "XX"
                else:
                    # get symbol by atomic number
                    symbol = GetAtomicSymbol(obatom.GetAtomicNum())
                try:
                    if read_fragments:
                        chain = str(obmol_id)
                    else:
                        chain = obatom.GetResidue().GetChain().strip()
                    atom = Atom(symbol, position, i, use_mass, chain)
                except (AttributeError, KeyError):
                    # molecule doesn't have chains, stuck with empty chains
                    atom = Atom(symbol, position, i, use_mass)
                adjacent = []
                for neighbour_atom in OBAtomAtomIter(obatom):
                    adjacent.append(neighbour_atom.GetIdx() - 1)
                atom.adjacent = MoleculeReader._remove_multi_bonds(adjacent)
                atoms.append(atom)

        if not len(atoms): #a stupid workaround for openbabel 3 bad handling of pdbs
            return None

        mol = Molecule(atoms)
        mol.metadata.file_content = mol_contents
        mol.metadata._title = title_contents
        return mol

    @staticmethod
    def read_xyz_connectivity(mol, conn_file):
        i = 0
        try:
            with open(conn_file, 'r') as file:
                for raw_line in file:
                    line = raw_line.split()
                    try:
                        atom_num = int(line.pop(0))
                    except (ValueError, IndexError):
                        raise ValueError("Input Error: Failed reading connectivity for atom " + str(i + 1))
                    if atom_num != i + 1:
                        raise ValueError("Input Error: Failed reading connectivity for atom " + str(i + 1))

                    neighbours = []
                    for neighbour_str in line:
                        try:
                            neighbour = int(neighbour_str) - 1  # Indexes in csm file start with 1
                        except ValueError:
                            raise ValueError("Input Error: Failed reading input for atom " + str(i + 1))
                        if neighbour >= len(mol):
                            raise ValueError("Input Error: Failed reading input for atom " + str(i + 1))
                        neighbours.append(neighbour)
                    mol.atoms[i].adjacent = MoleculeReader._remove_multi_bonds(neighbours)
                    i += 1
        except FileNotFoundError:
            raise FileNotFoundError("Failed to find connectivity file: " + str(conn_file))
        mol._create_bondset()

    @staticmethod
    def _read_csm_file(filename, ignore_sym=False, use_mass=False):
        """
        :param filename: Name of CSM file
        :param ignore_symbol: When true, the atom's symbol is not read
        :param use_mass: Use the atom's mass
        :return: A list of Atoms
        """

        with open(filename, 'r') as f:
            try:
                size = int(f.readline())
            except ValueError:
                raise ValueError("Input Error: Number of atoms not supplied")

            if size > 0:
                atoms = []
            else:
                return None

            for i in range(size):
                line = f.readline().split()
                try:
                    if ignore_sym:
                        symbol = "XX"
                    else:
                        symbol = line[0]
                    position = (float(line[1]), float(line[2]), float(line[3]))
                    atom = Atom(symbol, position, i, use_mass)
                except (ValueError, IndexError):
                    raise ValueError("Input Error: Failed reading input for atom " + str(i + 1))
                atoms.append(atom)

            for i in range(size):
                line = f.readline().split()
                try:
                    atom_num = int(line.pop(0))
                except (ValueError, IndexError):
                    raise ValueError("Input Error: Failed reading connectivity for atom " + str(i + 1))
                if atom_num != i + 1:
                    raise ValueError("Input Error: Failed reading connectivity for atom " + str(i + 1))

                neighbours = []
                for neighbour_str in line:
                    try:
                        neighbour = int(neighbour_str) - 1  # Indexes in csm file start with 1
                    except ValueError:
                        raise ValueError("Input Error: Failed reading input for atom " + str(i + 1))
                    if neighbour >= size:
                        raise ValueError("Input Error: Failed reading input for atom " + str(i + 1))
                    neighbours.append(neighbour)

                atoms[i].adjacent = MoleculeReader._remove_multi_bonds(neighbours)

            #chains = Chains() #used to get chain indices???
            try:
                numchains = int(f.readline())
                for i in range(numchains):
                    line = f.readline().split()
                    chain_name = line[0]
                    # chains[chain_name]=[]
                    for j in range(1, len(line)):
                        atom_index = int(line[j]) - 1
                        # chains[chain_name].append(atom_index)
                        atoms[atom_index]._chain = chains.index_map[chain_name]

            except:
                pass  # assume there are no chains
        mol = Molecule(atoms)
        with open(filename, 'r') as file:
            mol.metadata.file_content = file.read()
        return mol

    @staticmethod
    def _read_pdb_connectivity_and_chains(filename, mol, read_fragments, babel_bond, use_backbone=False):
        with open(filename, 'r') as file:
            # Count ATOM and HETATM lines, mapping them to our ATOM numbers.
            # In most PDBs ATOM 1 is our atom 0, and ATOM n is our n-1. However, in some cases
            # there are TER lines in the PDB, ATOMs after the TER are found at n-(TERCOUNT+1) in our list.
            # Some pdb files have multiple models even though we're only reading the first
            atom_map = {}
            cur_atom = 0
            # chains = Chains()


            for line in file:
                pdb_dict = PDBLine._pdb_line_to_dict(line)
                index = pdb_dict.atom_serial_number

                record_name = pdb_dict.record_name
                if record_name in ['ATOM', 'HETATM'] and cur_atom < len(mol):
                    if record_name == 'HETATM':
                        breakpt = 1
                    atom_map[index] = cur_atom
                    if use_backbone:
                        mol._atoms[cur_atom].atom_name = pdb_dict.atom_name  # change to the full symbol for use-backbone

                    chain_designation = pdb_dict.chain_id
                    if chain_designation != " ":
                        mol._atoms[cur_atom]._chain = str(chain_designation)

                    cur_atom += 1

                if not babel_bond:
                    if record_name == "CONECT":
                        # CONECT records are described here: http://www.bmsc.washington.edu/CrystaLinks/man/pdb/part_69.html
                        # After CONECT appears a series of 5 character atom numbers. There are no separating spaces in case
                        # the atom numbers have five digits, so we need to split the atom numbers differently

                        try:
                            fake_atom_index = pdb_dict.atom_serial_number
                            try:
                                atom_index = atom_map[fake_atom_index]
                            except KeyError:
                                if fake_atom_index > len(mol):
                                    continue
                            atom = mol._atoms[atom_index]

                            # add adjacency
                            adjacent = []
                            for adjacent_atom_index in pdb_dict.adjacent_atoms:
                                adjacent.append(atom_map[adjacent_atom_index])
                            atom.adjacent = MoleculeReader._remove_multi_bonds(adjacent)
                        except Exception as e:
                            raise ValueError("There was a problem reading connectivity from the pdb file." + str(e))
                    mol._create_bondset()  # refresh the bondset
        return mol

    @staticmethod
    def _create_pdb_with_sequence(mol, in_file_name, initialize=True,
                                  use_chains=False, babel_bond=False, read_fragments=False,
                                  ignore_hy=False, remove_hy=False, ignore_sym=False, use_mass=False,
                                  select_atoms=[], ignore_atoms=[], use_backbone=False):
        def read_atom(line, likeness_dict, cur_atom):
            pdb_dict = PDBLine._pdb_line_to_dict(line)
            atom_type = pdb_dict.atom_symbol
            remoteness = pdb_dict.remoteness
            sequence_number = pdb_dict.sequence_number
            key = tuple([atom_type, remoteness, sequence_number])
            if key not in likeness_dict:
                likeness_dict[key] = [cur_atom]
            else:
                likeness_dict[key].append(cur_atom)

        def set_equivalence_classes(mol, likeness_dict):
            groups = []
            for key in likeness_dict:
                groups.append(likeness_dict[key])

            mol._equivalence_classes = groups
            try:
                for group in groups:
                    for atom_index in group:
                        for equiv_index in group:
                            mol._atoms[atom_index].add_equivalence(equiv_index)
            except Exception as e:  # TODO: comment why this except is here (I don't actually remember)
                print(e)

        mol = MoleculeReader._read_pdb_connectivity_and_chains(in_file_name, mol, read_fragments, babel_bond, use_backbone)
        # mol.strip_atoms(remove_hy, select_atoms=select_atoms, ignore_atoms=ignore_atoms, use_backbone=use_backbone) # wait for answer from inbal, about the flags: select_atoms, ignore_atoms
        mol.strip_atoms(remove_hy, use_backbone=use_backbone)
        likeness_dict = {}
        cur_atom = 0

        with open(in_file_name, 'r') as file:
            for line in file:
                pdb_dict = PDBLine(line)
                if pdb_dict.record_name in ['ATOM', 'HETATM'] and cur_atom < len(mol):
                    if remove_hy or ignore_hy:
                        if pdb_dict.atom_symbol == "H":
                            continue
                    elif use_backbone and pdb_dict.atom_name not in ['C', 'CA', 'N', 'O']:
                        # pass when --use-backbone == True && the current atom is in the list of the backbone atoms
                        continue
                    read_atom(line, likeness_dict, cur_atom)
                    cur_atom += 1

        mol._create_bondset()
        mol.create_Q()

        # we have our own equivalence class function and hence cannot call the main initialization
        set_equivalence_classes(mol, likeness_dict)

        mol._initialize_chains(use_chains)
        if use_chains and len(set([len(mol.chains[i]) for i in mol.chains])) > 1:  # check if all the chains have the same length.
            raise ValueError("Error, all the chains of the molecule must have the same length, your chains' lengths: {}".format([len(mol.chains[i]) for i in mol.chains]))
        mol.normalize()

        return mol

    @staticmethod
    def redo_molecule(in_mol, **kwargs):
        '''
        when a comfile line argument includes a command that modifes a molecule, we resubmit molecules here to be redone
        :param in_mol:
        :param kwargs:
        :return:
        '''
        format = in_mol.metadata.format
        try:
            if kwargs["in_format"]:
                format = kwargs["in_format"]
        except KeyError:
            pass

        # question 1: do we need openbabel:
        # babel_bond
        try:
            babel_bond = kwargs["babel_bond"]
        except KeyError:
            babel_bond = in_mol.metadata.babel_bond

        obms = MoleculeReader._obm_from_strings(in_mol.metadata.file_content,
                                                format,
                                                babel_bond)
        if "read_fragments" in kwargs:
            out_mol = MoleculeReader.mol_from_obm(obms, format, **kwargs)
        else:
            obm = obms[0]
            out_mol = MoleculeReader.mol_from_obm([obm], format, **kwargs)

        kwargs.pop("in_file_name")

        out_mol = MoleculeReader._process_single_molecule(out_mol, in_mol.metadata.filepath, format, **kwargs)
        out_mol.metadata = in_mol.metadata
        out_mol.metadata.babel_bond = babel_bond
        return out_mol


def select_mols(mols, kwargs):
    select_mols=kwargs.get('select_mols', [])
    try:
        if select_mols:
            mols = [mols[i] for i in select_mols]
    except IndexError:
        raise IndexError("You have selected more molecules than you have input")
    return mols, select_mols


def mol_string_from_obm(obmol, format):
    obconversion = OBConversion()
    formatok = obconversion.SetOutFormat(format)
    if not formatok:
        raise ValueError("%s is not a recognised Open Babel format" %
                         format)
    return obconversion.WriteString(obmol)
