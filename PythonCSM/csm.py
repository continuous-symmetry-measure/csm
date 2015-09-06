import math

__author__ = 'zmbq'

"""
Performs some tests on the CSM C++ wrapper
"""
import sys
from input_output.writers import print_all_output
from calculations.preprocess_molecule import preprocess_molecule
from calculations.process_results import process_results
from calculations.csm_calculations_data import CSMCalculationsData
from calculations.csm_calculations import perform_operation
from arguments import process_arguments, create_parser
from CPP_wrapper import csm

MINDOUBLE = 1e-8
MAXDOUBLE = 100000000.0
APPROX_RUN_PER_SEC = 8e4


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
            if not 'perm' in csm_args and not 'dir' in csm_args:
                total_perms = csm.TotalNumberOfPemrutations()
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

        if csm_args['displayPerms']:
            perms = csm.GetMoleculePermutations()
            for i, perm in enumerate(perms):
                print("%5d: %s" % (i, perm))
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
                # TODO: options.opOrder = 2;
                data.opOrder = 2
                result = perform_operation(csm_args, data)

                if result.csm > MINDOUBLE:
                    data.operationType = "SN"
                    for i in range(2, csm_args['sn_max'] + 1, 2):
                        # TODO: options.opOrder = i
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
                # TODO: options.opOrder = chMinOrder;
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
        except:
            pass


if __name__ == '__main__':
    results = run_csm(sys.argv[1:])
