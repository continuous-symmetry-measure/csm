"""
Performs some tests on the CSM C++ wrapper
"""
import sys
from input_output.writers import print_all_output
from calculations.preprocess_molecule import preprocess_molecule
from calculations.process_results import process_results
from calculations.csm_calculations_data import CSMCalculationsData

__author__ = 'zmbq'

from arguments import process_arguments, create_parser
from CPP_wrapper import csm


def run_csm(args, print_output=True):
    try:
        parser = create_parser()
        result = parser.parse_args(args)
        csm_args = process_arguments(result)

        preprocess_molecule(csm_args)

        csm.SetCSMOptions(csm_args)

        # total_perms = csm.TotalNumberOfPemrutations()
        # print("Total number of permutations: %f" % total_perms)

        data = CSMCalculationsData(csm_args)
        perm_res = csm.RunSinglePerm(data)
        # TODO

        #results = csm.Calculate()

        results = {}
        results['molecule'] = perm_res.molecule
        results['outAtoms'] = perm_res.outAtoms
        results['dir'] = perm_res.dir
        results['csm'] = perm_res.csm
        results['dMin'] = perm_res.dMin
        results['perm'] = perm_res.perm
        results['localCSM'] = perm_res.localCSM
        results['operationType'] = perm_res.operationType

        process_results(results, csm_args)

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

