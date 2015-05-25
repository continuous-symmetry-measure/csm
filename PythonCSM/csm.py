"""
Performs some tests on the CSM C++ wrapper
"""
from input_output.writers import print_output

__author__ = 'zmbq'

from arguments import process_arguments, create_parser
from CPP_wrapper import csm
import sys

def process_input(args=None):
    if not args:
        args = sys.argv[1:]
    parser = create_parser()
    result = parser.parse_args(args)

    processed_args = process_arguments(result)

    if processed_args["ignoreHy"] or processed_args["removeHy"]:
        processed_args["obmol"].DeleteHydrogens()

    return processed_args

def run():
    try:
        args = process_input()
    except Exception as ex:
        print("Error parsing arguments: %s" % ex, file=sys.stderr)
        raise

    try:
        results = csm.RunCSM(args)
        # print_output(results, args)
    except Exception as ex:
        print("Error running CSM: %s" % ex, file=sys.stderr)
        if(args['writeOpenu']):
            print("ERR* %s *ERR" % ex)
        raise

if __name__=='__main__':
    print("PythonCSM is running")
    run()
    print("PythonCSM is done")

