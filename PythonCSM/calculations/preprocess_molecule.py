__author__ = 'YAEL'

from input_output.writers import print_equivalence_classes, print_molecule

# debug_file = open("debug_file.txt", "w")

def preprocess_molecule(csm_args):
    """
    Removes Hydrogen if needed and calculates equivalence classes
    :param csm_args:
    """

    # print_molecule(csm_args['molecule'], debug_file)

    if not csm_args['removeHy']:
        csm_args['equivalence_classes'] = find_equivalence_classes(csm_args['molecule'])
        # debug_file.write("Similarities after initSimilarity: there are %d groups\n" % len(csm_args['equivalence_classes']))
        # print_equivalence_classes(csm_args['equivalence_classes'], debug_file)

    if csm_args['ignoreHy'] or csm_args['removeHy']:
        csm_args["obmol"].DeleteHydrogens()
        remove_list = ["H", " H"]
        strip_atoms(csm_args, remove_list)
        # debug_file.write("Similarities after ignoreHy: there are %d groups\n" % len(csm_args['equivalence_classes']))
        # print_equivalence_classes(csm_args['equivalence_classes'],debug_file)

    # debug_file.close()


def find_equivalence_classes(atoms):
    group_num = 0
    groups = []
    atoms_size = len(atoms)
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
            if j in marked or len(atoms[i].adjacent) != len(atoms[j].adjacent) or atoms[i].symbol != atoms[j].symbol:
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

        for i,group in enumerate(groups):
            first_elem = group[0]
            group_size = len(group)

            sub_group = []

            # for each item in the group (except the first) check if it can be split or not
            for j in range(1,group_size):
                if not is_similar(atoms_group_num, atoms, group[j], first_elem):
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
        # print_equivalence_classes(groups, debug_file)

    #TODO: LOG(debug) << "Broken into groups with " << num_iters << " iterations.";
    print("Broken into groups with %d iterations." % num_iters)

    return groups


def is_similar(atoms_group_num, atoms, a, b):
    found = True
    mark = set()

    valency_a = len(atoms[a].adjacent)
    valency_b = len(atoms[b].adjacent)

    # for each of i's neighbours
    for i in range(valency_a):
        found = False

        for j in range(valency_b):
            if j in mark:
                continue

            if atoms_group_num[atoms[a].adjacent[i]] == atoms_group_num[atoms[b].adjacent[j]]:
                # the i-th neighbour of 'a' belongs to the same group as the j-th neighbour of 'b'
                found = True
                mark.add(j)
                break

        if not found:
            break

    return found


def strip_atoms(csm_args, remove_list):
    """
    Creates a new Molecule from m by removing atoms who's symbol is in the remove list
    :param csm_args:
    :param removeList: atomic symbols to remove
    """

    # find atoms in removeList
    to_remove = []
    size = len(csm_args['molecule'])
    for i in range(size):
        hits = 0
        for s in remove_list:
            if csm_args['molecule'][i].symbol == s:
                hits += 1
                break
        if hits > 0:
            to_remove.append(i)

    if len(to_remove) > 0:
        remove_atoms(csm_args, to_remove)


def remove_atoms(csm_args, to_remove):
    """
    Removes atoms with indexes in the to_remove list from the molecule
    :param csm_args:
    :param to_remove:
    """
    move_indexes = {}
    size = len(csm_args['molecule'])
    j = 0

    for i in range(size):
        if i == to_remove[j]:
            j += 1
        else:
            move_indexes[i] = i-j
    j -= 1

    # for i in move_indexes:
    #     debug_file.write("%d -> %d\n" % (i, move_indexes[i]))

    for i in range(size-1, 0, -1):
        if i == to_remove[j]:
            # remove the atom i
            csm_args['molecule'].pop(i)
            j -= 1
        else:
            # update the i-th atom adjacents
            l = len(csm_args['molecule'][i].adjacent)
            for k in range(l - 1, 0, -1):
                if csm_args['molecule'][i].adjacent[k] in move_indexes:
                    csm_args['molecule'][i].adjacent[k] = move_indexes[csm_args['molecule'][i].adjacent[k]]
                else:
                    csm_args['molecule'][i].adjacent.pop(k)

    if csm_args['ignoreHy']:
        # update indexes in equivalence classes
        groups_num = len(csm_args['equivalence_classes'])
        for i in range(groups_num - 1, -1, -1):
            group_size = len(csm_args['equivalence_classes'][i])
            for j in range(group_size - 1, -1, -1):
                if csm_args['equivalence_classes'][i][j] in move_indexes:
                    csm_args['equivalence_classes'][i][j] = move_indexes[csm_args['equivalence_classes'][i][j]]
                else:
                    csm_args['equivalence_classes'][i].pop(j)
            if len(csm_args['equivalence_classes'][i]) == 0:
                csm_args['equivalence_classes'].pop(i)

    if csm_args['removeHy']:
        csm_args['equivalence_classes'] = find_equivalence_classes(csm_args['molecule'])

