import csv
import json
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

def run(args=[]):
    print("CSM version %s" % csm.__version__)
    if not args:
        args = sys.argv[1:]
    csv_file = None
    try:
        # Read inputs
        dictionary_args = get_split_arguments(args)
        dictionary_args['molecule'], dictionary_args['perm'], dictionary_args['dirs'] = read_inputs(**dictionary_args)

        # Outputing permutations
        if dictionary_args['perms_csv_name']:
            csv_file = open(dictionary_args['perms_csv_name'], 'w')
            perm_writer = csv.writer(csv_file, lineterminator='\n')
            perm_writer.writerow(['Permutation', 'Direction', 'CSM'])
            exact_calculations.csm_state_tracer_func = lambda state: perm_writer.writerow([[p + 1 for p in state.perm],
                                                                                           state.dir,
                                                                                           state.csm, ])

        # run actual calculation
        if dictionary_args['calc_type'] == 'approx':
            result = approx_calculation(**dictionary_args)
        elif dictionary_args['calc_type'] == 'just_perms':
            result = perm_count(**dictionary_args)
        elif dictionary_args['calc_type'] == 'trivial':
            result = trivial_calculation(**dictionary_args)
        else:
            result = exact_calculation(**dictionary_args)


        print_results(result, dictionary_args)
        return result


    except Exception as e:
        if dictionary_args['json_output']:
            json_dict={
                "Error":str(e)
            }
            with open(dictionary_args['out_file_name'], 'w', encoding='utf-8') as f:
                json.dump(json_dict, f)
        raise

    finally:
        if csv_file:
            csv_file.close()


def run_no_return(args=[]):
    run(args)


def run_many_times():
    for i in range(100):
        print(i)
        run(args=sys.argv[1:])

if __name__ == '__main__':
    #timer = timeit.Timer(lambda: run_many_times())
    timer = timeit.Timer(lambda: run(args=sys.argv[1:]))
    time = timer.timeit(number=1)
    print("Runtime:", time)
