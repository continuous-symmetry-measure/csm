"""
Performs some tests on the CSM C++ wrapper
"""
from input_output.writers import print_all_output

__author__ = 'zmbq'

from arguments import process_arguments, create_parser
from CPP_wrapper import csm


if __name__=='__main__':
    parser = create_parser()
    result = parser.parse_args()  # Parse sys.args
    args = process_arguments(result)


    if args["ignoreHy"] or args["removeHy"]:
        args["obmol"].DeleteHydrogens()

    results = csm.RunCSM(args)

    print_all_output(results, args)



