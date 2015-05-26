"""
Performs some tests on the CSM C++ wrapper
"""
from input_output.writers import print_all_output
import sys

__author__ = 'zmbq'

from arguments import process_arguments, create_parser
from CPP_wrapper import csm


def parse_input(args):
    global openU

    parser = create_parser()
    result = parser.parse_args(args)

    return process_arguments(result)

def run(input_args):
    processed_args = None
    try:
        processed_args = parse_input(input_args)

        if processed_args["ignoreHy"] or processed_args["removeHy"]:
            processed_args["obmol"].DeleteHydrogens()

        results = csm.RunCSM(processed_args)
    except Exception as ex:
        print("CSM failed: %s" % str(ex), file=sys.stderr)
        if processed_args and processed_args['writeOpenu']:
            print("ERR* %s ERR*" % str(ex))
        raise

    print_all_output(results, processed_args)

if __name__=='__main__':
    try:
        run(sys.argv[1:])
    except ValueError as ve:
        raise
