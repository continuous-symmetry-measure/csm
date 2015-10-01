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
        csm_args['molecule'].find_equivalence_classes()
        # debug_file.write("Similarities after initSimilarity: there are %d groups\n" % len(csm_args['equivalence_classes']))
        # print_equivalence_classes(csm_args['equivalence_classes'], debug_file)

    if csm_args['ignoreHy'] or csm_args['removeHy']:
        if "obmol" in csm_args:
            csm_args["obmol"].DeleteHydrogens()
        remove_list = ["H", " H"]
        csm_args['molecule'].strip_atoms(remove_list, csm_args['ignoreHy'])
        # debug_file.write("Similarities after ignoreHy: there are %d groups\n" % len(csm_args['equivalence_classes']))
        # print_equivalence_classes(csm_args['equivalence_classes'],debug_file)

    # debug_file.close()

    csm_args['molecule'].normalize(csm_args['keepCenter'])

    #if csm_args['printPermutations']:
        # Print the equivalency classes
    #    print("Equivalence classes:")
    #    for cls in csm_args['molecule']._equivalence_classes:
    #        print("%3d: %s" % (len(cls), cls))


