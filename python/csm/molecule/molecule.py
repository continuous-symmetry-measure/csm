from collections import OrderedDict

from openbabel import OBAtomAtomIter, OBConversion, OBMol
from csm.molecule.atom import Atom, GetAtomicSymbol
from csm.molecule.normalizations import normalize_coords, de_normalize_coords
import logging
import numpy as np

logger = logging.getLogger("csm")


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
        return self._Q

    @property
    def atoms(self):
        return self._atoms

    @property
    def bondset(self):
        return self._bondset

    @property
    def equivalence_classes(self):
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
        return self._chains

    @property
    def obmol(self):
        return self._obmol

    def __len__(self):
        return len(self._atoms)

    def has_bond(self, atom_i, atom_j):
        if (atom_i, atom_j) in self._bondset:
            return True
        return False

    def atom_cords(self):
        atoms = []
        for atom in self._atoms:
            atoms.append(atom.pos)
        return atoms

    def _create_bondset(self):
        for i in range(len(self._atoms)):
            for match in self._atoms[i].adjacent:
                self._bondset.add((i, match))

    def _find_equivalence_classes(self):
        group_num = 0
        groups = []
        atoms_size = len(self._atoms)
        marked = set()

        atoms_group_num = {}

        logger.debug("Breaking molecule into similarity groups")

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

        logger.debug("initial number of groups:", group_num)
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
                    if not self.is_similar(atoms_group_num, group[j], first_elem):
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
        TODO: Improve this comment
        within each equivalence class, labels by chain
        self.group_chains=[ {key:[]}]
        """

        if not use_chains:
            self._chains = {'A': list(range(len(self.atoms)))}  # Default - one chain of all the atoms

        i=0
        self.chainkeys= {}
        for key in self._chains:
            self.chainkeys[key]=i
            i+=1

        group_chains=[]
        for group in self.equivalence_classes:
            chaingroup={}

            for i, atom_index in enumerate(group):
                try: #get chain-- if no chain, default to A
                    if use_chains:
                        chain= self.atoms[atom_index].chain
                    else:
                        chain= 'A'
                except:
                    chain= 'A'

                try:#add index to array in dict
                    chaingroup[self.chainkeys[chain]].append(atom_index)
                except: #create array in dict
                    chaingroup[self.chainkeys[chain]]=[atom_index]
            group_chains.append(chaingroup)
        self.group_chains = group_chains




    def is_similar(self, atoms_group_num, a, b):
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

    def _calculate_equivalency(self, remove_hy=False, ignore_hy=False):
        """
        Preprocess a molecule based on the arguments passed to CSM
        :param remove_hy: True if hydrogen atoms should be removed
        :param ignore_hy: True when hydrogen atoms should be ignored when calculating the equivalence classes
        :param kwargs: Place holder for all other csm_args.
        You can call it by passing **csm_args
        """
        if not remove_hy:
            self._find_equivalence_classes()

        if ignore_hy or remove_hy:
            if self._obmol:
                self._obmol.DeleteHydrogens()
            remove_list = ["H", " H"]
            self.strip_atoms(remove_list, ignore_hy)
            self._find_equivalence_classes()

    def strip_atoms(self, remove_list, ignore_hy):
            """
            Creates a new Molecule from m by removing atoms who's symbol is in the remove list
            :param csm_args:
            :param removeList: atomic symbols to remove
            """

            # find atoms in removeList
            to_remove = []
            size = len(self._atoms)
            for i in range(size):
                hits = 0
                for s in remove_list:
                    if self._atoms[i].symbol == s:
                        hits += 1
                        break
                if hits > 0:
                    to_remove.append(i)
            print(hits, "molecules of hydrogen removed or ignored")
            if len(to_remove) > 0:
                self.remove_atoms(to_remove, ignore_hy)

    def remove_atoms(self, to_remove, ignore_hy):
        """
        Removes atoms with indexes in the to_remove list from the molecule
        :param csm_args:
        :param to_remove:
        """
        move_indexes = {}
        size = len(self._atoms)
        j = 0

        for i in range(size):
            if i == to_remove[j]:
                j += 1
            else:
                move_indexes[i] = i - j
        j -= 1

        for i in range(size - 1, 0, -1):
            if i == to_remove[j]:
                # remove the atom i
                self._atoms.pop(i)
                if not ignore_hy:
                    if self.obmol:
                        self.obmol.DeleteAtom(self.obmol.GetAtom(i + 1))
                j -= 1
            else:
                # update the i-th atom adjacents
                l = len(self._atoms[i].adjacent)
                for k in range(l - 1, 0, -1):
                    if self._atoms[i].adjacent[k] in move_indexes:
                        self._atoms[i].adjacent[k] = move_indexes[self._atoms[i].adjacent[k]]
                    else:
                        self._atoms[i].adjacent.pop(k)

        if ignore_hy:
            # update indexes in equivalence classes
            groups_num = len(self._equivalence_classes)
            for i in range(groups_num - 1, -1, -1):
                group_size = len(self._equivalence_classes[i])
                for j in range(group_size - 1, -1, -1):
                    if self._equivalence_classes[i][j] in move_indexes:
                        self._equivalence_classes[i][j] = move_indexes[self._equivalence_classes[i][j]]
                    else:
                        self._equivalence_classes[i].pop(j)
                if len(self._equivalence_classes[i]) == 0:
                    self._equivalence_classes.pop(i)
        else:  # removeHy
            self._find_equivalence_classes()


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

    def create_Q(self):
        self._Q= np.array([np.array(atom.pos) for atom in self.atoms])


    def de_normalize(self):
        coords = [atom.pos for atom in self._atoms]
        denorm_coords = de_normalize_coords(coords, self.norm_factor)

        size = len(self._atoms)
        for i in range(size):
            self._atoms[i].pos = denorm_coords[i]

    def _complete_initialization(self, remove_hy, ignore_hy, use_chains):
        """
        Finish creating the molecule after reading the raw data
        """
        def diagnostics():
            lengths={}
            for group in self._equivalence_classes:
                try:
                    lengths[len(group)]+=1
                except:
                    lengths[len(group)]=1
            for key in lengths:
                print("%d group%s of length %d" %(lengths[key], 's' if lengths[key] else '', key))

            if use_chains:
                for chain in self.chains:
                    print ("Chain %s of length %d" % (chain, len(self.chains[chain])))
            else:
                print("--useChains not specified. using one simulated chain of len %d" %len(self.chains['A']))

        print("Breaking molecule into similarity groups")
        self._calculate_equivalency(remove_hy, ignore_hy)
        self._process_chains(use_chains)
        print("Broken into "+str(len(self._equivalence_classes))+" groups")
        diagnostics()
        self.normalize()

    @staticmethod
    def dummy_molecule(size):
        atoms=[]
        for i in range(size):
            atom = Atom("XX", (0,0,0), False)
            atoms.append(atom)
        mol= Molecule(atoms)
        group = [i for i in range(size)]
        mol._equivalence_classes=[group]
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
                  remove_hy=False, ignore_symm=False, use_mass=False, *args, **kwargs):
        if format == "csm":
            mol = Molecule._read_csm_file(in_file_name, ignore_symm, use_mass)
        else:
            obm = Molecule._obm_from_file(in_file_name, format, babel_bond)
            mol = Molecule._from_obm(obm, ignore_symm, use_mass)
        if initialize:
            mol._complete_initialization(remove_hy, ignore_hy, use_chains)
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
        conv = OBConversion()
        mol = OBMol()
        if not format:
            format = conv.FormatFromExt(filename)
            if not format:
                raise ValueError("Error discovering format from filename " + filename)
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
        chains = OrderedDict()
        for i in range(num_atoms):
            obatom = obmol.GetAtom(i + 1)
            if ignore_symm:
                symbol = "XX"
            else:
                # get symbol by atomic number
                symbol = GetAtomicSymbol(obatom.GetAtomicNum())
            position = (obatom.GetX(), obatom.GetY(), obatom.GetZ())
            chain = obatom.GetResidue().GetChain()
            if chain not in chains:
                chains[chain]=[]
            chains[chain].append(i)
            atom = Atom(symbol, position, use_mass, chain)
            adjacent = []
            iter = OBAtomAtomIter(obatom)
            for neighbour_atom in iter:
                adjacent.append(neighbour_atom.GetIdx() - 1)
            atom.adjacent = adjacent
            atoms.append(atom)

        equal=True
        for chain in chains:
            test= len(atoms)%len(chains[chain])
            if test!=0: #check that length of chain is, at minimum, a divisor of number of atoms
                equal=False
        #if not equal:
        #    raise Exception("Molecule's chains not of expected length: % num of chains, % num of molecules", (len(chains), len(atoms)))
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
                    atom = Atom(symbol, position, use_mass)
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

                atoms[i].adjacent = neighbours

            try:
                numchains=int(f.readline())
                chains=OrderedDict()
                for i in range(numchains):
                    line = f.readline().split()
                    chain_name= line[0]
                    chains[chain_name]=[]
                    for j in range(1, len(line)):
                        atoms[int(line[j])-1]._chain=chain_name
                        chains[chain_name].append(int(line[j])-1)
            except:
                pass
        return Molecule(atoms=atoms, chains=chains)
