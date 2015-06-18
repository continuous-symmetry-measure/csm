"""
Performs some tests on the CSM C++ wrapper
"""
import sys
from input_output.writers import print_all_output
from calculations.find_equivalence_classes import find_equivalence_classes
import os

__author__ = 'zmbq'

from arguments import process_arguments, create_parser
from CPP_wrapper import csm

def run_csm(args, print_output=True):
    try:
        parser = create_parser()
        result = parser.parse_args(args)
        csm_args = process_arguments(result)

        csm_args['equivalence_classes'] = find_equivalence_classes(csm_args['molecule'])
        if csm_args["ignoreHy"] or csm_args["removeHy"]:
            csm_args["obmol"].DeleteHydrogens()

        results = csm.RunCSM(csm_args)
        if print_output:
            print_all_output(results, csm_args)
        return results
    finally:
        try:
            csm_args['outFile'].close()
        except:
            pass

if __name__=='__main__':
    results = run_csm(sys.argv[1:])

