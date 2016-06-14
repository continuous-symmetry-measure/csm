print("csm_run.py imported")

class Stopwatch:
    def __init__(self):
        import time
        self._start = time.time()
        self._prev= self._start

    def elapsed(self):
        import time
        return time.time() - self._start

    def gap(self):
        import time
        gap= time.time()-self._prev
        self._prev=time.time()
        return gap

    def report(self, msg):
        return

stopwatch = Stopwatch()

stopwatch.report("Stopwatch started")
import csv
stopwatch.report("Imported csv")

import logging
stopwatch.report("Imported logging")

import sys
stopwatch.report("Imported sys")

import timeit
stopwatch.report("Imported timeit")

from csm.input_output.arguments import get_split_arguments
stopwatch.report("Import csm.input_output_arguments.get_split_arguments")

from csm.calculations.csm_calculations import exact_calculation, perm_count
stopwatch.report('imported csm.calculations_csmcalculations.exact_calculations, perm_count')

from csm.calculations.approx_calculations import approx_calculation, trivial_calculation
stopwatch.report('Imported approx_calculations')

from csm.calculations import csm_calculations
stopwatch.report("Imported csm_calculations")

from csm.input_output.readers import read_inputs
stopwatch.report('Imported read_inputs")'
                 '')
from csm.input_output.writers import print_results
stopwatch.report("Imported print_results")

import csm
stopwatch.report("Imported csm")

APPROX_RUN_PER_SEC = 8e4
sys.setrecursionlimit(10000)
stopwatch.report("Set recursion limit")

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
    stopwatch.report('run starting up')
    print("CSM version %s" % csm.__version__)
    if not args:
        args = sys.argv[1:] 
    csv_file = None
    try:
        # Read inputs
        in_args, calc_args, out_args = get_split_arguments(args)
        calc_args['molecule'], calc_args['perm'], calc_args['dir'] = read_inputs(**in_args)
        stopwatch.report("Arguments parsed")

        # logging:
        init_logging(**out_args)
        stopwatch.report("Logging initialized")

        # Outputing permutations
        if out_args['perms_csv_name']:
            csv_file = open(out_args['perms_csv_name'], 'w')
            perm_writer = csv.writer(csv_file, lineterminator='\n')
            perm_writer.writerow(['Permutation', 'Direction', 'CSM'])
            csm_calculations.csm_state_tracer_func = lambda state: perm_writer.writerow([[p + 1 for p in state.perm],
                                                                                         state.dir,
                                                                                         state.csm, ])

        # run actual calculation
        if calc_args['calc_type']=='approx':
            result = approx_calculation(**calc_args)
        elif calc_args['calc_type']=='just_perms':
            result = perm_count(**calc_args)
        elif calc_args['calc_type']=='trivial':
            result=trivial_calculation(**calc_args)
        else:
            result = exact_calculation(**calc_args)

        print_results(result, in_args, calc_args, out_args)
        return result

    finally:
        stopwatch.report("Inside run's finally")
        if csv_file:
            csv_file.close()

def run_no_return(args=[]):
    run(args)

if __name__ == '__main__':
    stopwatch.report("__name__ is __main__")
    timer = timeit.Timer(lambda: run(args=sys.argv[1:]))
    stopwatch.report("timeit initialized")
    time = timer.timeit(number=1)
    stopwatch.report("Time measured")
    print("Runtime:", time)
