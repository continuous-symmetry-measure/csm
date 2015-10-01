import math
import sys
__author__ = 'zmbq'

"""
Performs some tests on the CSM C++ wrapper
"""
from input_output.writers import print_all_output
from calculations.preprocess_molecule import preprocess_molecule
from calculations.process_results import process_results
from calculations.csm_calculations_data import CSMCalculationsData
from calculations.csm_calculations import perform_operation, MAXDOUBLE, total_number_of_permutations
from arguments import process_arguments, create_parser
from CPP_wrapper import csm

MINDOUBLE = 1e-8
APPROX_RUN_PER_SEC = 8e4


sys.setrecursionlimit(10000)


def run_csm(args, print_output=True):
    try:
        parser = create_parser()
        result = parser.parse_args(args)
        csm_args = process_arguments(result)

        preprocess_molecule(csm_args)

        csm.SetCSMOptions(csm_args)

        if csm_args['babelTest']:
            return None

        if not csm_args['findPerm']:
            if 'perm' not in csm_args and 'dir' not in csm_args:
                total_perms = total_number_of_permutations(csm_args)
                time = 1.0 * total_perms / 3600 / APPROX_RUN_PER_SEC
                if math.isnan(time):
                    # time is NaN
                    time = MAXDOUBLE
                print("Going to enumerate over %d permutations" % total_perms)
                print("Entire run should take approx. %.2lf hours on a 2.0Ghz Computer" % time)
            else:
                print("Using 1 permutation")
                print("Run should be instantaneous")

            if csm_args['timeOnly']:
                return None

        data = CSMCalculationsData(csm_args)

        # Code from the old mainRot.cpp
        if 'perm' in csm_args:
            result = csm.RunSinglePerm(data)
        else:
            if csm_args['type'] != 'CH':
                result = perform_operation(csm_args, data)
            else:
                # chirality support
                data.operationType = data.chMinType = "CS"
                data.opOrder = 2
                result = perform_operation(csm_args, data)

                if result.csm > MINDOUBLE:
                    data.operationType = "SN"
                    for i in range(2, csm_args['sn_max'] + 1, 2):
                        data.opOrder = i
                        ch_result = perform_operation(csm_args, data)
                        if ch_result.csm < result.csm:
                            result = ch_result
                            result.chMinType = 'SN'
                            result.chMinOrder = ch_result.opOrder

                        if result.csm < MINDOUBLE:
                            break

        if csm_args['printLocal']:
            if csm_args['type'] == 'CH':
                data.opOrder = result.chMinOrder
            local_res = csm.ComputeLocalCSM(data)
            result.localCSM = local_res.localCSM

        process_results(result, csm_args)

        if print_output:
            print_all_output(result, csm_args)

        return result
    finally:
        try:
            csm_args['outFile'].close()
            csm_args['outPermFile'].close()
        except:
            pass


if __name__ == '__main__':
    results = run_csm(sys.argv[1:])
