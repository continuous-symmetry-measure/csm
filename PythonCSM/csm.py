"""
Performs some tests on the CSM C++ wrapper
"""
__author__ = 'zmbq'

from arguments import process_arguments, create_parser
from CPP_wrapper import csm

if __name__=='__main__':
    parser = create_parser()
    result = parser.parse_args() # Parse sys.args
    args = process_arguments(result)

    if args["ignoreHy"] or args["removeHy"]:
        args["obmol"].DeleteHydrogens()

    output_dict = csm.RunCSM(args)

    # TODO: print output_dict
    # print_output(output_dict, args)