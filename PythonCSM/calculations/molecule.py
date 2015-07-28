from openbabel import OBAtom, OBElementTable

__author__ = 'zmbq'

_tbl = OBElementTable()


def GetAtomicMass(symbol):
    atomicNum = _tbl.GetAtomicNum(symbol)
    atom = OBAtom()
    atom.SetAtomicNum(atomicNum)
    atom.SetIsotope(0)
    return atom.GetAtomicMass()


def GetAtomicSymbol(atomic_num):
    return _tbl.GetSymbol(atomic_num)


class Atom:
    def __init__(self, symbol, pos, useMass=True):
        self._symbol = symbol
        self.adjacent = []
        self._pos = pos
        if useMass and symbol != 'XX':
            self._mass = GetAtomicMass(symbol)
        else:
            self._mass = 1.0

    @property
    def mass(self):
        return self._mass

    @property
    def pos(self):
        return self._pos

    @property
    def symbol(self):
        return self._symbol

    def __str__(self):
        return "Symbol: %s\tPos: %s\tAdjacent: %s" % (self.symbol, self.pos, self.adjacent)


class Molecule:
    # A Molecule has atoms and equivalency classes
    def __init__(self, atoms, equivalence_classes=None):
        self._atoms = atoms
        if equivalence_classes:
            self._equivalence_classes = equivalence_classes
        else:
            self._equivalence_classes = []

    @property
    def atoms(self):
        return self._atoms

    @property
    def equivalence_classes(self):
        return self._equivalence_classes


    @property
    def norm_factor(self):
        # Normalization factor. Defaults to 1.0 if normalize wasn't called
        return self._norm_factor

    def normalize(self):
        """  Normalize the molecule  """
        pass

    def find_equivalence_classes(self):
        group_num = 0
        groups = []
        atoms_size = len(self._atoms)
        marked = set()

        atoms_group_num = {}

        # TODO: LOG(debug) << "Breaking molecule into similarity groups";
        print("Breaking molecule into similarity groups")

        # break into initial groups by symbol and valency
        for i in range(atoms_size):
            if i in marked:
                continue

            groups.append([])

            for j in range(atoms_size):
                if j in marked\
                        or len(self._atoms[i].adjacent) != len(self._atoms[j].adjacent)\
                        or self._atoms[i].symbol != self._atoms[j].symbol:
                    continue

                groups[group_num].append(j)
                atoms_group_num[j] = group_num
                marked.add(j)

            group_num += 1

        # Comment from C++:
        # iteratively refine the breakdown into groups
        # In a previous version we had 'depth' iterations, this version breaks into subgroups at an infinite depth -
        # as long as there's something to break, it is broken

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
                        atoms_group_num[j]

                if len(sub_group) > 0:
                    divided_group = True
                    new_groups.append(sub_group)
                    for el in sub_group:
                        atoms_group_num[el] = group_num
                    group_num += 1
                    # remove elements of sub_group from group
                    groups[i] = [el for el in group if el != -1]

            groups.extend(new_groups)

        #TODO: LOG(debug) << "Broken into groups with " << num_iters << " iterations.";
        print("Broken into groups with %d iterations." % num_iters)

        self._equivalence_classes = groups

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
                move_indexes[i] = i-j
        j -= 1

        for i in range(size-1, 0, -1):
            if i == to_remove[j]:
                # remove the atom i
                self._atoms.pop(i)
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
            self.find_equivalence_classes()



