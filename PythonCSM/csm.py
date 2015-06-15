"""
Performs some tests on the CSM C++ wrapper
"""
from input_output.writers import print_all_output
from calculations.find_equivalence_classes import find_equivalence_classes

__author__ = 'zmbq'

from arguments import process_arguments, create_parser
from CPP_wrapper import csm


if __name__=='__main__':
    parser = create_parser()
    result = parser.parse_args()  # Parse sys.args
    args = process_arguments(result)
    # args['equivalence_classes'] = csm.CallInitSimilarity(args['molecule'])
    args['equivalence_classes'] = find_equivalence_classes(args['molecule'])

    if args["ignoreHy"] or args["removeHy"]:
        args["obmol"].DeleteHydrogens()

    results = csm.RunCSM(args)

    print_all_output(results, args)

