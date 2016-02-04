import logging
import sys

from a_calculations.permuters import MoleculeLegalPermuter
from input_output.arguments import get_split_arguments
from a_calculations.csm_calculations_data import CSMCalculationsData
from a_calculations.csm_calculations import exact_calculation
from input_output.readers import read_inputs
from input_output.writers import print_results

APPROX_RUN_PER_SEC = 8e4
sys.setrecursionlimit(10000)

logger = None


def init_logging(log_file_name=None, *args, **kwargs):
    global logger

    if log_file_name:
        logging.basicConfig(filename=log_file_name, level=logging.DEBUG,
                            format='[%(asctime)-15s] [%(levelname)s] [%(name)s]: %(message)s')
    else:
        logging.basicConfig(level=logging.ERROR)
    logger = logging.getLogger("csm")


def run_csm(args={}):
    # Read inputs
    in_args, calc_args, out_args = get_split_arguments(args)
    calc_args['molecule'], calc_args['perm'], calc_args['dir'] = read_inputs(**in_args)
    calc_args['permuter_class'] = MoleculeLegalPermuter

    print("Running calculation with arguments: ", calc_args)
    # logging:
    init_logging(**out_args)

    # run actual calculation
    if calc_args['find_perm']:
        raise NotImplementedError("No approx yet")
#        result = approx_calculation(**calc_args)
    else:
        result = exact_calculation(**calc_args)

    print_results(result, in_args, calc_args, out_args)

if __name__ == '__main__':
    results = run_csm(args=sys.argv[1:])
