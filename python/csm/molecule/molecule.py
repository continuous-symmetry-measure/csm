"""
@author: Devora Witty
"""

from collections import OrderedDict
import copy
from openbabel import OBAtomAtomIter, OBConversion, OBMol, OBMolAtomIter
from csm.molecule.atom import Atom, GetAtomicSymbol
from csm.molecule.normalizations import normalize_coords, de_normalize_coords, calculate_norm_factor
import logging
import numpy as np

logger = logging.getLogger("csm")



class Chains(OrderedDict):
    '''
    two sets of keys, string keys and integer keys that match the index of the chain
    the inner dict is composed of integer keys, string keys must be translated to int before sending on to inner dict
    '''
    def __init__(self):
        self._indexes_to_strings=[]
        self._strings_to_indexes={}
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
            key=str(key)
        self._indexes_to_strings.append(key)
        index=self._indexes_to_strings.index(key)
        self._strings_to_indexes[key]=index
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
        zipped_arr=[(key, self.__getitem__(key)) for key in self._indexes_to_strings]
        return zipped_arr

    def from_array(self, arr_of_tuples):
        for key, val in arr_of_tuples:
            self.__setitem__(key, val)

    def index_to_string(self, index):
        return self._indexes_to_strings[index]

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
            self._norm_factor = 1.0 #is recalculated by molecule._complete_initialization

            self._calculate_center_of_mass()
            self._create_bondset()
            self.create_Q()

        # not included:
        #equivalence class initialization, chain initialization
            self._equivalence_classes = []
            self._chains = Chains()
            self._groups_with_internal_chains=[]
            self._chains_with_internal_groups={}
            self._chain_equivalences=[]

        #getting rid of obmol:
            #self._obmol = obmol
            self._filename=None
            self._deleted_atom_indices=[]
            self._format=None
            self._babel_bond=False


    def copy(self):
        #deepcopy is used only for atoms,
        # because atoms are the only property changed between runs of Directions that necessitated copying
        # (due to call to denormalize)

        m=Molecule(to_copy=True)
        m._atoms=copy.deepcopy(self.atoms)
        m._norm_factor=self.norm_factor

        #m._obmol=self.obmol
        m._filename=self._filename
        m._deleted_atom_indices=self._deleted_atom_indices
        m._format= self._format
        m._babel_bond = self._babel_bond

        m._bondset=self.bondset
        m._Q=self.Q
        m._center_of_mass=self._center_of_mass

        m._equivalence_classes=self.equivalence_classes
        m._chains=self._chains
        m._groups_with_internal_chains=self._groups_with_internal_chains
        m._chains_with_internal_groups=self._chains_with_internal_groups
        m._chain_equivalences=self._chain_equivalences
        return m

    def to_json(self):
        return {
            #critical to include
            "atoms":[atom.to_json() for atom in self.atoms],

            # trivial to include
            "norm factor": self.norm_factor,
            "center of mass":self.center_of_mass,

            #expensive to recalculate
            "equivalence classes": self.equivalence_classes, #the most expensive part of loading large molecule
            "groups_with_internal_chains":self.groups_with_internal_chains,
            "chains_with_internal_groups":self.chains_with_internal_groups,
            "chain_equivalences":self.chain_equivalences,

            #Classes:
            #obmol: needed for printing:
            "deleted indices":self._deleted_atom_indices,
            "filename":self._filename,
            "format":self._format,
            "babel_bond": self._babel_bond,
            #chains
            "chains": self.chains.to_array()

            # unsure whether worth the bother of including
            # "bondset":self.bondset,
            # "Q":self.Q, #almost definitely not worth the bother
                }

    @staticmethod
    def from_json(json):
        atoms=[Atom.from_json(a) for a in json["atoms"]]
        m=Molecule(atoms)

        m.norm_factor=json["norm factor"]

        c=Chains()
        c.from_array(json["chains"])
        m._chains=c

        m._filename=json["filename"]
        m._deleted_atom_indices = json["deleted indices"]
        m._format=json["format"]
        m._babel_bond=json["babel_bond"]

        m._center_of_mass=json["center of mass"]
        m._equivalence_classes=json["equivalence classes"]
        unfixed_groups=json["groups_with_internal_chains"]
        m._groups_with_internal_chains=[]
        for group in unfixed_groups:
            fixed={int(key): value for key, value in group.items()}
            m._groups_with_internal_chains.append(fixed)
        unfixed_chains=json["chains_with_internal_groups"]
        fixed = {int(key): value for key, value in unfixed_chains.items()}
        m._chains_with_internal_groups=fixed
        m._chain_equivalences=json["chain_equivalences"]

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
        self._bondset=set()
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

        #create chains, and overwrite atom.chain with index of chain rather than string
        chains=Chains()
        for i, atom in enumerate(self.atoms):
            chain=atom.chain
            if not use_chains:
                chain='Simulated chain'
            if chain not in chains:
                chains[chain] = []
            chains[chain].append(i)
            atom.chain=chains._strings_to_indexes[chain]
        self._chains=chains

        #within each equivalence class, labels by chain,
        # within each chain, group by equivalence class
        groups_with_internal_chains=[]
        num_equiv=len(self.equivalence_classes)
        chains_with_internal_groups={chain:[None] * num_equiv for chain in self.chains}

        for group_index, group in enumerate(self.equivalence_classes):
            chaingroup={}
            for atom_index in group:
                chain_index= self.atoms[atom_index].chain

                # add index to array in chaingroup dict
                try:
                    chaingroup[chain_index].append(atom_index)
                except KeyError:
                    chaingroup[chain_index]=[atom_index]
                #add index to classed_chains
                try:
                     chains_with_internal_groups[chain_index][group_index].append(atom_index)
                except AttributeError:
                    chains_with_internal_groups[chain_index][group_index]=[atom_index]
            #add the completed chaingroup dict to group_chains arrat
            groups_with_internal_chains.append(chaingroup)

        self._groups_with_internal_chains = groups_with_internal_chains
        self._chains_with_internal_groups = chains_with_internal_groups

        #Calculates chain equivalencies. (e.g, chains A and B are equivalent, chains C and D are equivalent)
        #we define chains as equivalent if the number of atoms they have in each equivalence class are the same
        chain_equivalences=[]
        marked=[False]*len(self.chains)

        for chain_index, chainkey in enumerate(self.chains):
            if marked[chain_index]:
                continue
            equiv=[]
            for chain2_index, chainkey2 in enumerate(self.chains):
                #start with a simple length check to spare checking equivalence classes if chains arent same length to begin with
                #if len(self.chains[chainkey])!=len(self.chains[chainkey2]):
                #    continue
                same_lengths=True
                for group in groups_with_internal_chains:
                    if not same_lengths:
                        break
                    try:
                        length = len(group[chain_index])
                    except KeyError: #that chain is notn in this equivalence group
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
        self._chain_equivalences=chain_equivalences



    def _calculate_equivalency(self, remove_hy, ignore_hy):
        """
        Preprocess a molecule based on the arguments passed to CSM
        :param remove_hy: True if hydrogen atoms should be removed
        :param ignore_hy: True when hydrogen atoms should be ignored when calculating the equivalence classes
        :param kwargs: Place holder for all other csm_args.
        You can call it by passing **csm_args
        """
        if ignore_hy or remove_hy:
            self.strip_atoms(remove_hy, ignore_hy)

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

        logger.debug("initial number of groups:"+ str(group_num))
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

        logger.debug("Broken into %d groups with %d iterations." % (group_num, num_iters))

        self._equivalence_classes = groups
        for group in groups:
            for atom_index in group:
                for equiv_index in group:
                    self._atoms[atom_index].add_equivalence(equiv_index)

    def strip_atoms(self, remove_hy, ignore_hy):
        """
            Creates a new Molecule from m by removing atoms who's symbol is in the remove list
            :param csm_args:
            :param removeList: atomic symbols to remove
        """

        remove_list = ["H", " H"]
        removed_atoms=[]
        fixed_indexes=[i for i in range(len(self))]

        for i in range(len(self._atoms)):
            if self._atoms[i].symbol in remove_list:
                removed_atoms.append(i)
                fixed_indexes[i]=None
            else:
                # however many atoms have been removed up to this index is the amount its index needs adjusting by
                fixed_indexes[i]-=len(removed_atoms)

        #adjust the indices before we do any popping whatsoever
        for i, atom in enumerate(self._atoms):
            adjacent_new=[]
            for adjacent in atom.adjacent:
                if fixed_indexes[adjacent] is not None:
                    adjacent_new.append(fixed_indexes[adjacent])
            atom.adjacent=adjacent_new

        for to_remove in reversed(removed_atoms): #reversed order because popping changes indexes after
            self._atoms.pop(to_remove)
            self._deleted_atom_indices.append(to_remove)
            #if remove_hy: #this is meant to affect print at end
                #self._obmol.DeleteAtom(self._obmol.GetAtom(to_remove + 1))

        logger.debug(len(removed_atoms), "molecules of hydrogen removed or ignored")

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

    def de_normalize(self):
        coords = [atom.pos for atom in self._atoms]
        denorm_coords = de_normalize_coords(coords, self.norm_factor)

        size = len(self._atoms)
        for i in range(size):
            self._atoms[i].pos = denorm_coords[i]
        self.create_Q()

    def create_Q(self):
        self._Q= np.array([np.array(atom.pos) for atom in self.atoms])

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
            print("%d group%s of length %d" % (lengths[key], 's' if lengths[key] and lengths[key] > 1 else '', key))

        if display_chains:
            if len(self.chains)>1:
                for chain in self.chains:
                    print("Chain %s of length %d" % (self.chains._indexes_to_strings[chain], len(self.chains[chain])))
                print("%d equivalence class%s of chains" % (len(self.chain_equivalences), 'es' if lengths[key] else ''))
                for chaingroup in self.chain_equivalences:
                    chainstring = "Group of length " + str(len(chaingroup)) + ":"
                    for index in chaingroup:
                        chainstring += " "
                        chainstring += str(self.chains._indexes_to_strings[index])
                    print(str(chainstring))
            else:
                print("Molecule has no chains")

    def _complete_initialization(self, use_chains, remove_hy, ignore_hy=False ):
        """
        Finish creating the molecule after reading the raw data
        """
        self._calculate_equivalency(remove_hy, ignore_hy)
        self._initialize_chains(use_chains)
        self.normalize()

    @staticmethod
    def xyz_string(atoms, positions=None, header=""):
        '''
        :param atoms: an array of atoms
        :param positions: optional, an array of positions to overwrite positions in the array of atoms
        :return:
        '''
        coords=str(len(atoms))
        coords+="\n"+header+"\n"
        for i, atom in enumerate(atoms):
            coords+=str(atom.symbol)
            for index, coor in enumerate(atom.pos):
                coords+="\t"
                if positions is not None:
                    coords+=str(positions[i][index])
                else:
                    coords+=str(coor)
            coords+="\n"
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
        atoms=[]

        for i in range(size):
            atom = Atom("C", (0,0,0), i, False)
            atoms.append(atom)
        mol= Molecule(atoms)
        mol._equivalence_classes=groups
        return mol

    @staticmethod
    def dummy_molecule_from_coords(coords, groups=None):
        atoms=[]

        if groups is None:
            groups=[[i for i in range(len(coords))]]

        for fake_atom, group in enumerate(groups):
            for atom_index in group:
                coord=coords[atom_index]
                atom=Atom(str(fake_atom), coord, atom_index, False)
                atoms.append(atom)

        mol= Molecule(atoms)
        mol._calculate_equivalency(False, False)
        mol.normalize()
        return mol

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
                    remove_hy=False, ignore_symm=False, use_mass=False):
        """
        create a Molecule from a string
        :param string: the string to create the molecule from
        :param format: the format of the string (any BabelBond supported format, eg "mol", "xyz")
        :param initialize: boolean, default True, when True equivalence classes and chains are calculated for the molecule
        :param use_chains: boolean, default False, when True chains are read from the string 
        :param babel_bond: boolean, default False, when True OpenBabel will attempt to guess connectivity information
        :param remove_hy: boolean, default False, when True hydrogen atoms will be removed
        :param ignore_symm: boolean, default False, when True all atom symbols will be "X" instead of elements 
                            (this affects equivalence class calculation)
        :param use_mass: boolean, default False, when True atomic mass will be used, when False all atoms have a mass of 1
        :return: an instance of class Molecule
        """
        # note: useMass is used when creating molecule, even though it is actually about creating the normalization
        # step one: get the molecule object
        obm = MoleculeReader._obm_from_string(string, format, babel_bond)
        mol = MoleculeReader._from_obm(obm, ignore_symm, use_mass)


        if initialize:
            mol._complete_initialization(use_chains, remove_hy)

        return mol

    @staticmethod
    def from_file(in_file_name, format=None, initialize=True,
                  use_chains=False, babel_bond=False,
                  remove_hy=False, ignore_symm=False, use_mass=False,
                  read_fragments=False, use_sequence=False,
                  keep_structure=False,
                  *args, **kwargs):
        """
        :param in_file_name: the name of the file to read the molecule from
        :param format: the format of the string (any BabelBond supported format, eg "mol", "xyz")
        :param initialize: boolean, default True, when True equivalence classes and chains are calculated for the molecule
        :param use_chains: boolean, default False, when True chains are read from the string 
        :param babel_bond: boolean, default False, when True OpenBabel will attempt to guess connectivity information
        :param remove_hy: boolean, default False, when True hydrogen atoms will be removed
        :param ignore_symm: boolean, default False, when True all atom symbols will be "X" instead of elements 
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
        def get_format(form):
            if format.lower()=="csm":
                return "csm"
            conv = OBConversion()
            if not form:
                form = conv.FormatFromExt(in_file_name)
                if not form:
                    raise ValueError("Error discovering format from filename " + in_file_name)
            return form

        def set_obmol_field(mol):
            mol._filename=in_file_name
            mol._babel_bond=babel_bond
            mol._format=format

        format = get_format(format)

        if use_sequence:
            if format.lower() != 'pdb':
                raise ValueError("--use-sequence only works with PDB files")

            mol = MoleculeReader._create_pdb_with_sequence(in_file_name, format, use_chains=use_chains, babel_bond=babel_bond,
                            read_fragments=read_fragments, remove_hy=remove_hy, ignore_symm=ignore_symm, use_mass=use_mass)
            #we initialize mol from within pdb_with_sequence because otherwise equivalnce classes would be overwritten
            set_obmol_field(mol)
            return mol

        if format == "csm":
                mol = MoleculeReader._read_csm_file(in_file_name, ignore_symm, use_mass)

        else:
                obm = MoleculeReader._obm_from_file(in_file_name, format, babel_bond)
                mol = MoleculeReader._from_obm(obm, ignore_symm, use_mass, read_fragments)
                if format=="pdb":
                    mol=MoleculeReader._read_pdb_connectivity_and_chains(in_file_name, mol, read_fragments, babel_bond)
                if not mol.bondset:
                    if keep_structure:
                        raise ValueError("User input --keep-structure but input molecule has no bonds. Did you forget --babel-bond?")
                    else:
                        logger.warn("Input molecule has no bond information")

        if initialize:
            mol._complete_initialization(use_chains, remove_hy)
            if len(mol.chains)<2:
                if read_fragments:
                    logger.warn("Although you input --read-fragments, no fragments were found in file. "
                            "Fragments are marked by $$$ in mol files or by model/endmdl in pdb files")
                elif use_chains:
                    logger.warn("You specified --use-chains but molecule only has one chain")
        set_obmol_field(mol)
        return mol


    @staticmethod
    def _obm_from_string(string, format, babel_bond=None):
        conv = OBConversion()
        obmol = OBMol()
        if not conv.SetInFormat(format):
            raise ValueError("Error setting openbabel format to" + format)
        if not babel_bond:
            conv.SetOptions("b", conv.INOPTIONS)
        conv.ReadString(obmol, string)
        return [obmol]

    @staticmethod
    def _obm_from_file(filename, format=None, babel_bond=None):
        """
        :param filename: name of file to open
        :param format: molecule format of file (eg xyz, pdb)
        :param babelBond:
        :return:
        """
        obmol = OBMol()
        conv = OBConversion()
        if not conv.SetInFormat(format):
            raise ValueError("Error setting openbabel format to" + format)
        if not babel_bond:
            conv.SetOptions("b", conv.INOPTIONS)
        notatend = conv.ReadFile(obmol, filename)
        if not notatend:
            raise ValueError("Error reading file " + filename + " using OpenBabel")
        obmols=[]
        while notatend:
            obmols.append(obmol)
            obmol = OBMol()
            notatend = conv.Read(obmol)

        return obmols

    @staticmethod
    def _from_obm(obmols, ignore_symm=False, use_mass=False, read_fragments=False):
        """
        :param obmol: OBmol molecule
        :param args_dict: dictionary of processed command line arguments
        :return: A list of Atoms and a list of chains
        """
        if not read_fragments:
            obmols=obmols[:1]
        atoms = []
        for obmol_id, obmol in enumerate(obmols):
            for i, obatom in enumerate(OBMolAtomIter(obmol)):
                position = (obatom.GetX(), obatom.GetY(), obatom.GetZ())
                if ignore_symm:
                    symbol = "XX"
                else:
                    # get symbol by atomic number
                    symbol = GetAtomicSymbol(obatom.GetAtomicNum())
                try:
                    if read_fragments:
                        chain=str(obmol_id)
                    else:
                        chain = obatom.GetResidue().GetChain().strip()
                    atom = Atom(symbol, position, i, use_mass, chain)
                except (AttributeError, KeyError):
                    #molecule doesn't have chains, stuck with empty chains
                    atom = Atom(symbol, position, i, use_mass)
                adjacent = []
                for neighbour_atom in OBAtomAtomIter(obatom):
                    adjacent.append(neighbour_atom.GetIdx() - 1)
                atom.adjacent =  MoleculeReader._remove_multi_bonds(adjacent)
                atoms.append(atom)

        mol = Molecule(atoms)
        return mol

    @staticmethod
    def _read_csm_file(filename, ignore_symbol=False, use_mass=False):
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
                    if ignore_symbol:
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

                atoms[i].adjacent =  MoleculeReader._remove_multi_bonds(neighbours)

            #chains = Chains() #used to get chain indices???
            try:
                numchains=int(f.readline())
                for i in range(numchains):
                    line = f.readline().split()
                    chain_name= line[0]
            #        chains[chain_name]=[]
                    for j in range(1, len(line)):
                        atom_index=int(line[j]) - 1
            #            chains[chain_name].append(atom_index)
                        atoms[atom_index]._chain=chains.index_map[chain_name]

            except:
                pass #assume there are no chains
        return Molecule(atoms)


    @staticmethod
    def _read_pdb_connectivity_and_chains(filename, mol, read_fragments, babel_bond):
        with open(filename, 'r') as file:
            # Count ATOM and HETATM lines, mapping them to our ATOM numbers.
            # In most PDBs ATOM 1 is our atom 0, and ATOM n is our n-1. However, in some cases
            # there are TER lines in the PDB, ATOMs after the TER are found at n-(TERCOUNT+1) in our list.
            #Some pdb files have multiple models even though we're only reading the first
            atom_map = {}
            cur_atom = 0
            #chains = Chains()

            for line in file:
                try:
                    index = int(line[6:11])
                except (ValueError, IndexError):
                    index = None

                if line[0:6] in ['ATOM  ','HETATM'] and cur_atom<len(mol):
                    if line[0:6]=='HETATM':
                        breakpt=1
                    atom_map[index] = cur_atom

                    #handle chains
                    #ATOMxxxxxx1xxNxxxLYSxA
                    #HETATMxxxx3xxHxxxHOHxxxxx1
                    if line[0:6]=='ATOM  ':
                        chain_designation=line[21]
                    if line[0:6]=='HETATM':
                        if line[21]!=" ":
                            chain_designation = line[21]
                        else:
                            chain_designation = line[25]
                    if chain_designation!=" ":
                        mol._atoms[cur_atom]._chain = str(chain_designation)

                    cur_atom += 1

                if not babel_bond:
                    if line[0:6] == "CONECT":
                        # CONECT records are described here: http://www.bmsc.washington.edu/CrystaLinks/man/pdb/part_69.html
                        # After CONECT appears a series of 5 character atom numbers. There are no separating spaces in case
                        # the atom numbers have five digits, so we need to split the atom numbers differently

                        try:
                            line = line.strip()  # Remove trailing whitespace
                            fake_atom_index = int(line[6:11])
                            try:
                                atom_index = atom_map[fake_atom_index]
                            except KeyError:
                                if fake_atom_index>len(mol):
                                    continue
                            atom = mol._atoms[atom_index]


                            #add adjacency
                            adjacent = []
                            for i in range(11, len(line), 5):
                                adjacent_atom_index = int(line[i:i + 5])
                                adjacent.append(atom_map[adjacent_atom_index])
                            atom.adjacent = MoleculeReader._remove_multi_bonds(adjacent)
                        except Exception as e:
                            raise ValueError("There was a problem reading connectivity from the pdb file." + str(e))
                    mol._create_bondset() #refresh the bondset
        return mol


    @staticmethod
    def _create_pdb_with_sequence(in_file_name, format="pdb", initialize=True,
                                  use_chains=False, babel_bond=False, read_fragments=False,
                                  ignore_hy=False, remove_hy=False, ignore_symm=False, use_mass=False):
        def read_atom(line, likeness_dict, cur_atom):
            atom_type = line[12:14]
            remoteness = line[14]
            serial_number = line[22:26]
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

        obm = MoleculeReader._obm_from_file(in_file_name, format, babel_bond)
        mol = MoleculeReader._from_obm(obm, ignore_symm, use_mass)
        mol = MoleculeReader._read_pdb_connectivity_and_chains(in_file_name, mol, read_fragments, babel_bond)
        if remove_hy or ignore_hy:
            mol.strip_atoms(remove_hy, ignore_hy)
        likeness_dict = {}
        cur_atom=0


        with open(in_file_name, 'r') as file:
            for line in file:
                if line[0:6] in ['ATOM  ','HETATM'] and cur_atom<len(mol):
                    if remove_hy or ignore_hy:
                        if line[12:14] in ['H', ' H', 'H ']:
                            continue
                    read_atom(line, likeness_dict, cur_atom)
                    cur_atom+=1

        mol._create_bondset()
        mol.create_Q()

        #we have our own equivalence class function and hence cannot call the main initialization
        set_equivalence_classes(mol, likeness_dict)

        mol._initialize_chains(use_chains)
        mol.normalize()

        return mol

