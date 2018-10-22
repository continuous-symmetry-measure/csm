"""
@author: Devora Witty
"""
import sys
import os
from collections import OrderedDict
import copy
from openbabel import OBAtomAtomIter, OBConversion, OBMol, OBMolAtomIter, obErrorLog, obError
from csm.molecule.atom import Atom, GetAtomicSymbol
from csm.molecule.normalizations import normalize_coords, de_normalize_coords, calculate_norm_factor
import logging
import numpy as np
import json
from csm.input_output.formatters import csm_log as print

logger = logging.getLogger("csm")

ob_debug = False
if not ob_debug:
    obErrorLog.SetOutputLevel(obError)


def get_format(format, filename):
    if not format:
        format = filename.split(".")[-1]
    if format.lower() == "csm":
        return "csm"
    conv = OBConversion()
    if not conv.SetInFormat(format):
        format = conv.FormatFromExt(filename)
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

    def index_to_string(self, index):
        return self._indexes_to_strings[index]

class MoleculeMetaData:
    '''
    This class is primarily used to store metadata needed to write results, although format+filecontent+babel_bond
    are sometimes used to recreate a molecule from scratch, see: redo_molecule

    format, filename, and babel_bond are set in _initialize_single_molecule
    file_content is set in read_obm or read_csm
    index is set in read_multiple_molecules, or, revoltingly, in do_commands after calling redo_molecule
    '''
    def __init__(self, mol_contents=[], format=None, filename="", babel_bond=False, index=0, title=None):
        self.file_content = mol_contents
        self.format = format
        self.filename=filename
        self.babel_bond = babel_bond
        self.index=index
        self.title=title
        self.use_filename=True

    @staticmethod
    def from_dict(self, dict):
        file_content = dict["file_content"]
        format = dict["format"]
        filename=dict["filename"]
        babel_bond = dict["babel_bond"]
        index=dict["index"]
        title=dict["title"]
        m= MoleculeMetaData(file_content, format, filename, babel_bond, index, title)
        m.use_filename=dict["use_filename"]

    def to_dict(self):
        return {
            "file_content":self.file_content,
            "format":self.format,
            "filename":self.filename,
            "babel_bond":self.babel_bond,
            "index":self.index,
            "title":self.title,
            "use_filename":self.use_filename
        }

    def header(self, no_file_format=False, no_str_format=False):
        if self.use_filename:
            if no_file_format:
                return self.filename[:-4]
            return self.filename

        mol_index=self.index+1 #start from 1 instead of 0

        if self.title: #replace with internal index if relevant
            start_index=self.title.find("mol_index=")
            if start_index != -1:
                end_index=self.title.find(";")
                start_index=start_index+10
                mol_index=int(self.title[start_index:end_index])
        if no_str_format:
            return str(mol_index)

        mol_str="%04d" % mol_index
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
            self.has_been_normalized=None

            self.metadata=MoleculeMetaData()


    def copy(self):
        # deepcopy is used only for atoms,
        # because atoms are the only property changed between runs of Directions that necessitated copying
        # (due to call to denormalize)

        m = Molecule(to_copy=True)
        m._atoms = copy.deepcopy(self.atoms)
        m._norm_factor = self.norm_factor

        m._deleted_atom_indices = self._deleted_atom_indices
        m.metadata=self.metadata

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
            "metadata":self.metadata.to_dict(),

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
        m.has_been_normalized=in_dict["normalized"]

        c = Chains()
        c.from_array(in_dict["chains"])
        m._chains = c

        m._deleted_atom_indices = in_dict["deleted indices"]
        m.metadata= MoleculeMetaData.from_dict(in_dict["metadata"])

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
        :param remove_hy: True if hydrogen atoms should be removed
        :param ignore_hy: True when hydrogen atoms should be ignored when calculating the equivalence classes
        :param kwargs: Place holder for all other csm_args.
        You can call it by passing **csm_args
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

        #logger.debug("initial number of groups:" + str(group_num))
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

        #logger.debug("Broken into %d groups with %d iterations." % (group_num, num_iters))

        self._equivalence_classes = groups
        for group in groups:
            for atom_index in group:
                for equiv_index in group:
                    self._atoms[atom_index].add_equivalence(equiv_index)

    def strip_atoms(self, remove_hy=False, select_atoms=[]):
        """
            Creates a new Molecule from m by removing atoms who's symbol is in the remove list
            :param csm_args:
            :param removeList: atomic symbols to remove
        """

        removed_atoms = []
        fixed_indexes = [i for i in range(len(self))]

        for i in range(len(self._atoms)):
            if remove_hy:
                if self._atoms[i].symbol =="H":
                    removed_atoms.append(i)
                    fixed_indexes[i] = None
                else:
                    # however many atoms have been removed up to this index is the amount its index needs adjusting by
                    fixed_indexes[i] -= len(removed_atoms)

            if select_atoms:
                if i not in select_atoms:
                    removed_atoms.append(i)
                    fixed_indexes[i] = None
                else:
                    # however many atoms have been removed up to this index is the amount its index needs adjusting by
                    fixed_indexes[i] -= len(removed_atoms)


        # adjust the indices before we do any popping whatsoever
        for i, atom in enumerate(self._atoms):
            adjacent_new = []
            for adjacent in atom.adjacent:
                if fixed_indexes[adjacent] is not None:
                    adjacent_new.append(fixed_indexes[adjacent])
            atom.adjacent = adjacent_new

        for to_remove in reversed(removed_atoms):  # reversed order because popping changes indexes after
            self._atoms.pop(to_remove)
            self._deleted_atom_indices.append(to_remove)
            # if remove_hy: #this is meant to affect print at end
            # self._obmol.DeleteAtom(self._obmol.GetAtom(to_remove + 1))

        #logger.debug(len(removed_atoms), "atoms removed")

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
        self.has_been_normalized=True

    def de_normalize(self):
        coords = [atom.pos for atom in self._atoms]
        denorm_coords = de_normalize_coords(coords, self.norm_factor)

        size = len(self._atoms)
        for i in range(size):
            self._atoms[i].pos = denorm_coords[i]
        self.create_Q()
        self.has_been_normalized=False

    def create_Q(self):
        self._Q = np.array([np.array(atom.pos) for atom in self.atoms])

    def print_equivalence_class_summary(self, display_chains=True):
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
            print("%d group%s of length %d" % (lengths[key], 's' if lengths[key] and lengths[key] > 1 else '', key))

        if True: #len(self.chains) > 1:
                for chain in self.chains:
                    print("Chain %s of length %d" % (self.chains._indexes_to_strings[chain], len(self.chains[chain])))
                print("%d equivalence class%s of chains" % (len(self.chain_equivalences), 'es' if lengths[key] else ''))
                for chaingroup in self.chain_equivalences:
                    chainstring = "Group of length " + str(len(chaingroup)) + ":"
                    for index in chaingroup:
                        chainstring += " "
                        chainstring += str(self.chains._indexes_to_strings[index])
                    print(str(chainstring))
            #else:
            #    print("Molecule has no chains")

    def _complete_initialization(self, use_chains, remove_hy, select_atoms=[]):
        """
        Finish creating the molecule after reading the raw data
        """
        self.strip_atoms(remove_hy, select_atoms)
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
        mol._calculate_equivalency(False, False)
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
                self.remoteness=self.atom_symbol[1]
                self.branch_designation=self.atom_symbol[2]
                self.duplicate_index=self.atom_symbol[3]
            except: #TODO: check that what information is available gets recorded
                pass #because it may be attempting to access information that doesn't exist so it's ok...
        else:
            self.atom_symbol=pdb_line[12:14].strip()
            self.remoteness = pdb_line[14]
            self.branch_designation = pdb_line[15]
            self.alternate_location = pdb_line[16]

        if not self.atom_symbol:
            raise ValueError("PDB missing atom symbol information in atom name (columns 13-16). Please check that your pdb is formatted correctly")

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
    def from_string(string, format, initialize=True,
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
            mol = MoleculeReader.mol_from_obm(obm, format, ignore_sym=ignore_sym, use_mass=use_mass, read_fragments=read_fragments)
        return MoleculeReader._process_single_molecule(mol, in_file_name, format, initialize,
                                                       use_chains, babel_bond,
                                                       remove_hy, ignore_sym, use_mass,
                                                       read_fragments, use_sequence,
                                                       keep_structure, select_atoms, conn_file)

    @staticmethod
    def _process_single_molecule(mol, in_file_name, format, initialize=True,
                                 use_chains=False, babel_bond=False,
                                 remove_hy=False, ignore_sym=False, use_mass=False,
                                 read_fragments=False, use_sequence=False,
                                 keep_structure=False, select_atoms=[], conn_file=None, **kwargs):

        mol.metadata.format=format
        mol.metadata.babel_bond=babel_bond
        mol.metadata.filename=os.path.basename(in_file_name)

        if use_sequence:
            if format.lower() == 'pdb':
                mol = MoleculeReader._create_pdb_with_sequence (mol, in_file_name, use_chains=use_chains,
                                                           babel_bond=babel_bond,
                                                           read_fragments=read_fragments, remove_hy=remove_hy,
                                                           ignore_sym=ignore_sym, use_mass=use_mass)
                # we initialize mol from within pdb_with_sequence because otherwise equivalnce classes would be overwritten
                return mol

        if format == "pdb":
            mol = MoleculeReader._read_pdb_connectivity_and_chains(in_file_name, mol, read_fragments, babel_bond)
        if conn_file and format=="xyz":
            MoleculeReader.read_xyz_connectivity(mol, conn_file)
        if not mol.bondset:
            if keep_structure:
                raise ValueError(
                    "User input --keep-structure but input molecule has no bonds. Did you forget --babel-bond?")
            else:
                print("Warning: Input molecule has no bond information")

        if initialize:
            mol._complete_initialization(use_chains, remove_hy, select_atoms)
            if len(mol.chains) < 2:
                raise ValueError("The molecule you have provided has less than 2 chains. "
                                 "Protein-csm only works with molecules with 2 or more chains.")
        return mol

    @staticmethod
    def multiple_from_file(in_file_name, in_format=None, initialize=True,
                           use_chains=False, babel_bond=False,
                           remove_hy=False, ignore_sym=False, use_mass=False,
                           read_fragments=False, use_sequence=False,
                           keep_structure=False, select_atoms=[], conn_file=None,
                           *args, **kwargs):
        mols = []
        format = get_format(in_format, in_file_name)

        if format == "csm":
            mol = MoleculeReader._read_csm_file(in_file_name, ignore_sym, use_mass)
            mols.append(mol)

        else:
            obms = MoleculeReader._obm_from_file(in_file_name, format, babel_bond)
            if read_fragments:
                mol = MoleculeReader.mol_from_obm(obms, format, babel_bond=babel_bond, ignore_sym=ignore_sym, use_mass=use_mass, read_fragments=read_fragments)
                mols.append(mol)

            else:
                for obm in obms:
                    mol = MoleculeReader.mol_from_obm([obm], format, babel_bond=babel_bond, ignore_sym=ignore_sym, use_mass=use_mass)
                    mols.append(mol)

        processed_mols=[]
        use_filename=True
        if len(mols)>1:
            use_filename=False
        for index, mol in enumerate(mols):
            p_mol = MoleculeReader._process_single_molecule(mol, in_file_name, format, initialize,
                                                          use_chains, babel_bond,
                                                          remove_hy, ignore_sym, use_mass,
                                                          read_fragments, use_sequence,
                                                          keep_structure, select_atoms, conn_file)
            p_mol.metadata.index=index
            p_mol.metadata.use_filename=use_filename
            processed_mols.append(p_mol)
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
        obmols=[]
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
        if format=="txt":
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
    def mol_from_obm(obmols, format, babel_bond=False, ignore_sym=False, use_mass=False, read_fragments=False, **kwargs):
        """
        :param obmol: OBmol molecule
        :param args_dict: dictionary of processed command line arguments
        :return: A list of Atoms and a list of chains
        """
        if not read_fragments:
            obmols = obmols[:1]
        atoms = []
        mol_contents=[]
        for obmol_id, obmol in enumerate(obmols):
            title_contents=obmol.GetTitle()
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

        mol = Molecule(atoms)
        mol.metadata.file_content=mol_contents
        mol.metadata.title=title_contents
        return mol

    @staticmethod
    def _read_pdb_connectivity_and_chains(filename, mol, read_fragments, babel_bond):
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
                                  ignore_hy=False, remove_hy=False, ignore_sym=False, use_mass=False):
        def read_atom(line, likeness_dict, cur_atom):
            pdb_dict = PDBLine._pdb_line_to_dict(line)
            atom_type = pdb_dict.atom_symbol
            remoteness = pdb_dict.remoteness
            serial_number = pdb_dict.sequence_number
            key = tuple([atom_type, remoteness, serial_number])
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

        mol = MoleculeReader._read_pdb_connectivity_and_chains(in_file_name, mol, read_fragments, babel_bond)
        if remove_hy or ignore_hy:
            mol.strip_atoms(remove_hy, ignore_hy)
        likeness_dict = {}
        cur_atom = 0

        with open(in_file_name, 'r') as file:
            for line in file:
                pdb_dict = PDBLine(line)
                if pdb_dict.record_name in ['ATOM', 'HETATM'] and cur_atom < len(mol):
                    if remove_hy or ignore_hy:
                        if pdb_dict.atom_symbol=="H":
                            continue
                    read_atom(line, likeness_dict, cur_atom)
                    cur_atom += 1

        mol._create_bondset()
        mol.create_Q()

        # we have our own equivalence class function and hence cannot call the main initialization
        set_equivalence_classes(mol, likeness_dict)

        mol._initialize_chains(use_chains)
        mol.normalize()

        return mol




def mol_string_from_obm(obmol, format):
    obconversion = OBConversion()
    formatok = obconversion.SetOutFormat(format)
    if not formatok:
        raise ValueError("%s is not a recognised Open Babel format" %
                         format)
    return obconversion.WriteString(obmol)