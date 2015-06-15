__author__ = 'YAEL'
import numpy as np

def find_equivalence_classes(atoms):
    group_num = 0
    groups = []
    atoms_size = len(atoms)
    marked = set()

    atoms_group_num = {}
    # atoms_group_num = np.zeros(atoms_size, dtype=int)

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

    #TODO: LOG(debug) << "Broken into groups with " << num_iters << " iterations.";
    print("Broken into groups with %d iterations." % num_iters)

    return groups


def is_similar(atoms_group_num, atoms, a, b):
    found = True
    atoms_size = len(atoms)
    mark = set() # [False for i in range(atoms_size)]

    valency_a = len(atoms[a].adjacent)
    valency_b = len(atoms[b].adjacent)

    # for each of i's neighbours
    for i in range(valency_a):
        found = False

        for j in range(valency_b):
            #if mark[j]:
            if j in mark:
                continue

            if atoms_group_num[atoms[a].adjacent[i]] == atoms_group_num[atoms[b].adjacent[j]]:
                # the i-th neighbour of 'a' belongs to the same group as the j-th neighbour of 'b'
                found = True
                mark.add(j)
                # mark[j] = True
                break

        if not found:
            break

    return found
