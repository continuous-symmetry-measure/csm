from collections import OrderedDict

from openbabel import OBAtomAtomIter, OBConversion, OBMol
from csm.molecule.atom import Atom, GetAtomicSymbol
from csm.molecule.normalizations import normalize_coords, de_normalize_coords
import logging
import numpy as np

logger = logging.getLogger("csm")


def remove_multi_bonds(bonds):
    l= list(set(bonds))
    l.sort()
    return l

class Chains(OrderedDict):
    def __init__(self):
        self.str_keys=[]
        self.index_map={}
        super().__init__()
    def __getitem__(self, key):
        if isinstance(key, str):
            index_key=self.index_map[key]
            return super().__getitem__(index_key)
        if isinstance(key, int):
            return super().__getitem__(key)

        raise ValueError
    def __setitem__(self, key, value):
        if isinstance(key, str):
            self.str_keys.append(key)
            index_key=self.str_keys.index(key)
            self.index_map[key]=index_key
            return super().__setitem__(index_key, value)
        raise ValueError

    def __contains__(self, item):
        if isinstance(item, int):
            return item in self.__dict__
        if isinstance(item, str):
            try:
                index_key=self.index_map[item]
                return super().__contains__(index_key)
            except:
                return False





class Molecule:
    def __init__(self, atoms={}, chains={}, norm_factor=1.0, obmol=None):

        self._atoms = atoms
        self._chains = chains
        self._bondset = set()
        self._equivalence_classes = []
        self._norm_factor = norm_factor
        self._flags = {}
        self._create_bondset()
        self._obmol = obmol
        self.create_Q()

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

    def set_norm_factor(self, nf):
        self._norm_factor = nf

    @property
    def chains(self):
        '''
        :return: a dictionary of the molecule's chains, with keys being the name of the chains
        '''
        return self._chains

    @property
    def size(self):
        return len(self._atoms)

    @property
    def obmol(self):
        return self._obmol

    def __len__(self):
        return len(self._atoms)

    def _create_bondset(self):
        for i in range(len(self._atoms)):
            for match in self._atoms[i].adjacent:
                self._bondset.add((i, match))

    def _find_equivalence_classes(self):
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



    def _process_chains(self, use_chains):
        """
        1. if there are no chains, creates one large simulated chain
        3. within each equivalence class, labels by chain
            self.group_chains=[array of equivalence classes:
                            {dictionary of chains:
                                [array of indexes belonging to that chain in that equivalence class]
                            }
                        ]
        4. Calculates chain equivalencies. (e.g, chains A and B are equivalent, chains C and D are equivalent)
        """

        #1. if there are no chains, creates one large simulated chain
        try:
            if not use_chains:
                self._chains=Chains()
                self._chains['Sim'] = list(range(len(self.atoms)))  # Default - one chain of all the atoms, "simulated"
                #print("--no-chains specified. Using one simulated chain of len %d" % len(self._chains['Sim']))
            else:
                check=self.atoms[0].chain #see if there even are chains
        except:
            self._chains = Chains()
            self._chains['Sim'] = list(range(len(self.atoms)))
            print("Molecule has no chains specified. Using one simulated chain of len %d" % len(self.chains[0]))



        #3. within each equivalence class, labels by chain,
        # within each chain, group by equivalence class
        groups_with_internal_chains=[]
        num_equiv=len(self.equivalence_classes)
        chains_with_internal_groups={chain:[None] * num_equiv for chain in self.chains}

        for group_index, group in enumerate(self.equivalence_classes):
            chaingroup={}
            for atom_index in group:
                # get chain-- if no chains, or not use_chains default to chain 0
                try:
                    if use_chains:
                        chain_index= self.atoms[atom_index].chain
                    else:
                        chain_index= 0
                except AttributeError:
                    chain_index= 0

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
        self.groups_with_internal_chains = groups_with_internal_chains
        self.chains_with_internal_groups = chains_with_internal_groups

        #4. Calculates chain equivalencies. (e.g, chains A and B are equivalent, chains C and D are equivalent)
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
                    except:
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
        self.chain_equivalences=chain_equivalences



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

        self._find_equivalence_classes()

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
            if remove_hy: #this is meant to affect print at end
                self._obmol.DeleteAtom(self._obmol.GetAtom(to_remove + 1))

        logger.debug(len(removed_atoms), "molecules of hydrogen removed or ignored")



    def normalize(self):
        """
        Normalize the molecule
        """
        coords = [atom.pos for atom in self._atoms]
        masses = [atom.mass for atom in self._atoms]
        (norm_coords, self._norm_factor) = normalize_coords(coords, masses)
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


    def display_info(self, display_chains):
        """
        Displays information about equivalence classes and chains
        """
        lengths = {}
        for group in self._equivalence_classes:
            try:
                lengths[len(group)] += 1
            except:
                lengths[len(group)] = 1
        for key in lengths:
            print("%d group%s of length %d" % (lengths[key], 's' if lengths[key] and lengths[key] > 1 else '', key))

        if display_chains:
            for chain in self.chains:
                print("Chain %s of length %d" % (self.chains.str_keys[chain], len(self.chains[chain])))
            print("%d equivalence class%s of chains" % (len(self.chain_equivalences), 'es' if lengths[key] else ''))
            for chaingroup in self.chain_equivalences:
                chainstring = "Group of length " + str(len(chaingroup)) + ":"
                for index in chaingroup:
                    chainstring += " "
                    chainstring += str(self.chains.str_keys[index])
                print(str(chainstring))

    def _complete_initialization(self, remove_hy, ignore_hy, use_chains):
        """
        Finish creating the molecule after reading the raw data
        """

        print("Breaking molecule into equivalence class groups")
        self._calculate_equivalency(remove_hy, ignore_hy)
        #print("Broken into " + str(len(self._equivalence_classes)) + " groups")
        self._process_chains(use_chains)
        try:
            self.display_info(use_chains)
        except:
            pass
        self.normalize()

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
    def molecule_from_coords(coords, groups=None):
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

    @staticmethod
    def from_string(string, format, initialize=True, use_chains=False, babel_bond=False, ignore_hy=False,
                    remove_hy=False, ignore_symm=False, use_mass=False):
        # note: useMass is used when creating molecule, even though it is actually about creating the normalization

        # step one: get the molecule object
        obm = Molecule._obm_from_string(string, format, babel_bond)
        mol = Molecule._from_obm(obm, ignore_symm, use_mass)
        if initialize:
            mol._complete_initialization(remove_hy, ignore_hy, use_chains)

        return mol

    @staticmethod
    def from_file(in_file_name, initialize=True, format=None, use_chains=False, babel_bond=False, ignore_hy=False,
                  remove_hy=False, ignore_symm=False, use_mass=False, keep_structure=False, no_babel=False, use_sequence=False,
                  *args, **kwargs):
        def get_format(form):
            if format.lower()=="csm":
                return "csm"
            conv = OBConversion()
            if not form:
                form = conv.FormatFromExt(in_file_name)
                if not form:
                    raise ValueError("Error discovering format from filename " + in_file_name)
            return form

        format = get_format(format)

        if use_sequence:
            if format.lower() != 'pdb':
                raise ValueError("--use-sequence only works with PDB files")
            mol = Molecule.create_pdb_with_sequence(in_file_name, initialize, format, use_chains, babel_bond, ignore_hy,
                                                    remove_hy, ignore_symm, use_mass, keep_structure, no_babel)
            return mol



        else:
            if format == "csm":
                mol = Molecule._read_csm_file(in_file_name, ignore_symm, use_mass)

            else:
                obm = Molecule._obm_from_file(in_file_name, format, babel_bond)
                mol = Molecule._from_obm(obm, ignore_symm, use_mass)
                if format=="pdb" and not babel_bond:
                    mol=Molecule._read_pdb_connectivity(in_file_name, mol)
                if not mol.bondset:
                    if keep_structure:
                        raise ValueError("User input --keep-structure but input molecule has no bonds. Did you forget --babel-bond?")
                    if babel_bond:
                        print("Molecule file does not have connectivity information. Using babelbond to create bonds.")
                        obm = Molecule._obm_from_file(in_file_name, format, True)
                        mol = Molecule._from_obm(obm, ignore_symm, use_mass)
                    else:
                        logger.warn("Input molecule has no bond information")

        if initialize:
            mol._complete_initialization(remove_hy, ignore_hy, use_chains)

        return mol

    @staticmethod
    def _read_pdb_connectivity(filename, mol):
        with open(filename, 'r') as file:
            # Count ATOM and HETATM lines, mapping them to our ATOM numbers.
            # In most PDBs ATOM 1 is our atom 0, and ATOM n is our n-1. However, in some cases
            # there are TER lines in the PDB, ATOMs after the TER are found at n-2 in our list.
            atom_map = {}
            chain_map={}
            cur_atom = 0
            chains=Chains()

            for line in file:
                parts = line.split()
                try:
                    index = int(parts[1])
                except (ValueError, IndexError):
                    index = None

                if parts[0] in ['ATOM', 'HETATM']:
                    atom_map[index] = cur_atom
                    chain_map[index] = line[25:26]
                    cur_atom += 1

                if line[0:6] == "CONECT":
                    # CONECT records are described here: http://www.bmsc.washington.edu/CrystaLinks/man/pdb/part_69.html
                    # After CONECT appears a series of 5 character atom numbers. There are no separating spaces in case
                    # the atom numbers have five digits, so we need to split the atom numbers differently
                    if False:
                        line = line.strip()  # Remove trailing whitespace
                        first_atom_index = int(line[6:11])
                        first_atom = mol._atoms[atom_map[first_atom_index]]

                        adjacent = []
                        for i in range(11, len(line), 5):
                            adjacent_atom_index = int(line[i:i+5])
                            adjacent.append(atom_map[adjacent_atom_index])
                        first_atom.adjacent = remove_multi_bonds(adjacent)
                    try:
                        line = line.strip()  # Remove trailing whitespace
                        fake_atom_index = int(line[6:11])
                        atom_index=atom_map[fake_atom_index]
                        atom = mol._atoms[atom_index]

                        #add chains
                        chain=chain_map[fake_atom_index]
                        if chain not in chains:
                            chains[chain] = []
                        chains[chain].append(atom_index)
                        mol._atoms[atom_index]._chain=chains.index_map[chain] #we store atom chain as int
                    except Exception as e:
                        raise ValueError("There was a problem reading connectivity from the pdb file." + str(e))
        mol._chains=chains
        return mol

    @staticmethod
    def create_pdb_with_sequence(filename, initialize=True, format=None, use_chains=False, babel_bond=False, ignore_hy=False,
                                 remove_hy=False, ignore_symm=False, use_mass=False, keep_structure=False, no_babel=False):

        def read_atom(line, likeness_dict, index):
            record_name=line[0:6]
            if record_name in ["ATOM  ", "HETATM"]:
                #handle equivalence class:
                atom_type=line[12:14]
                remoteness=line[14]
                serial_number=line[22:26]
                key=tuple([atom_type, remoteness, serial_number])
                val=index
                index+=1
                if key not in likeness_dict:
                    likeness_dict[key]=[val]
                else:
                    likeness_dict[key].append(val)

                #handle chains:
                chain=line[21]
            return index

        def set_equivalence_classes(mol, likeness_dict):
            groups=[]
            for key in likeness_dict:
                groups.append(likeness_dict[key])

            mol._equivalence_classes = groups
            try:
                for group in groups:
                    for atom_index in group:
                        for equiv_index in group:
                            mol._atoms[atom_index].add_equivalence(equiv_index)
            except:
                pass

        print("Breaking molecule into equivalency class groups based on protein sequence")
        obm = Molecule._obm_from_file(filename, format, babel_bond)
        mol = Molecule._from_obm(obm, ignore_symm, use_mass)

        likeness_dict={}
        index=0

        with open(filename, 'r') as file:
            for line in file:
                index=read_atom(line, likeness_dict, index)

        mol._create_bondset()
        set_equivalence_classes(mol, likeness_dict)
        print(len(likeness_dict), "equivalence groups")
        mol._process_chains(use_chains)
        mol.display_info(use_chains)
        mol.normalize()

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
        return obmol

    @staticmethod
    def _obm_from_file(filename, format=None, babel_bond=None):
        """
        :param filename: name of file to open
        :param format: molecule format of file (eg xyz, pdb)
        :param babelBond:
        :return:
        """
        mol = OBMol()
        conv = OBConversion()
        if not conv.SetInFormat(format):
            raise ValueError("Error setting openbabel format to" + format)
        if not babel_bond:
            conv.SetOptions("b", conv.INOPTIONS)
        if not conv.ReadFile(mol, filename):
            raise ValueError("Error reading file " + filename + " using OpenBabel")
        return mol

    @staticmethod
    def _from_obm(obmol, ignore_symm=False, use_mass=False):
        """
        :param obmol: OBmol molecule
        :param args_dict: dictionary of processed command line arguments
        :return: A list of Atoms and a list of chains
        """
        num_atoms = obmol.NumAtoms()
        atoms = []
        chains = Chains()
        for i in range(num_atoms):
            obatom = obmol.GetAtom(i + 1)
            if ignore_symm:
                symbol = "XX"
            else:
                # get symbol by atomic number
                symbol = GetAtomicSymbol(obatom.GetAtomicNum())
            position = (obatom.GetX(), obatom.GetY(), obatom.GetZ())
            try:
                chain = obatom.GetResidue().GetChain().strip()
                if chain:
                    if chain not in chains:
                        chains[chain]=[]
                chains[chain].append(i)
                atom = Atom(symbol, position, i, use_mass, chains.index_map[chain])
            except (AttributeError, KeyError):
                #molecule doesn't have chains, stuck with empty chains
                atom = Atom(symbol, position, i, use_mass)
            adjacent = []
            iter = OBAtomAtomIter(obatom)
            for neighbour_atom in iter:
                adjacent.append(neighbour_atom.GetIdx() - 1)
            atom.adjacent =  remove_multi_bonds(adjacent)
            atoms.append(atom)

        equal=True
        for chain in chains:
            test= len(atoms)%len(chains[chain])
            if test!=0: #check that length of chain is, at minimum, a divisor of number of atoms
                equal=False

        mol = Molecule(atoms=atoms, chains=chains, obmol=obmol)
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

                atoms[i].adjacent =  remove_multi_bonds(neighbours)

            chains = Chains()
            try:
                numchains=int(f.readline())
                for i in range(numchains):
                    line = f.readline().split()
                    chain_name= line[0]
                    chains[chain_name]=[]
                    for j in range(1, len(line)):
                        atom_index=int(line[j]) - 1
                        chains[chain_name].append(atom_index)
                        atoms[atom_index]._chain=chains.index_map[chain_name]

            except:
                pass #assume there are no chains
        return Molecule(atoms=atoms, chains=chains)


    def print_debug(self):
        '''
        matches same function in C++, used for comparing printouts
        '''
        #print("Equivalent Groups:")
        #for item in self.equivalence_classes:
        #    print(item)

        print("========DEBUG INFORMATION========")

        #print("molecule:")
        #print("size=", self.size)
        #for atom in self.atoms:
        #    print(atom.symbol, atom.pos)

        print("\nconnectivity:")
        for i in range(self.size):
            mystr = str(i+1)
            for adj in self.atoms[i].adjacent:
                mystr+=" "
                mystr+=str(adj+1)
            mystr+="\t\tTotal="
            mystr+=str(len(self.atoms[i].adjacent))
            print(mystr)

        #print("\nComplete Equivalent Groups for each atom:")
        #for i in range(self.size):
        #    mystr = str(i+1)
        #    for adj in self.atoms[i].equivalency:
        #        mystr+=" "
        #        mystr+=str(adj+1)
        #    print(mystr)

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


