import csv
import logging
import sys
import timeit
from csm.input_output.arguments import get_split_arguments
from csm.calculations.exact_calculations import exact_calculation, perm_count
from csm.calculations.approx_calculations import approx_calculation, trivial_calculation
from csm.calculations import exact_calculations
from csm.input_output.readers import read_inputs
from csm.input_output.writers import print_results
import csm

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


def run(args=[]):
    print("CSM version %s" % csm.__version__)
    if not args:
        args = sys.argv[1:]
    csv_file = None
    try:
        # Read inputs
        in_args, calc_args, out_args = get_split_arguments(args)
        calc_args['molecule'], calc_args['perm'], calc_args['dirs'] = read_inputs(**in_args)

        # logging:
        init_logging(**out_args)

        # Outputing permutations
        if out_args['perms_csv_name']:
            csv_file = open(out_args['perms_csv_name'], 'w')
            perm_writer = csv.writer(csv_file, lineterminator='\n')
            perm_writer.writerow(['Permutation', 'Direction', 'CSM'])
            exact_calculations.csm_state_tracer_func = lambda state: perm_writer.writerow([[p + 1 for p in state.perm],
                                                                                           state.dir,
                                                                                           state.csm, ])

        # run actual calculation
        if calc_args['calc_type'] == 'approx':
            result = approx_calculation(**calc_args)
        elif calc_args['calc_type'] == 'just_perms':
            result = perm_count(**calc_args)
        elif calc_args['calc_type'] == 'trivial':
            result = trivial_calculation(**calc_args)
        else:
            result = exact_calculation(**calc_args)


        print_results(result, in_args, calc_args, out_args)
        return result

    except:
        raise

    finally:
        if csv_file:
            csv_file.close()


def run_no_return(args=[]):
    run(args)


if __name__ == '__main__':
    timer = timeit.Timer(lambda: run(args=sys.argv[1:]))
    time = timer.timeit(number=1)
    print("Runtime:", time)
