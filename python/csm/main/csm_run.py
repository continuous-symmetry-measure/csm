import csv
import json
import logging
import sys
import timeit

from csm.input_output.arguments import get_split_arguments
from csm.calculations import Approx, Trivial, Exact
from csm.input_output.readers import read_inputs
from csm.input_output.writers import print_results, FileWriter
from csm import __version__
from csm.molecule.molecule import MoleculeReader

sys.setrecursionlimit(10000)


def run(args=[]):
    print("CSM version %s" % __version__)
    if not args:
        args = sys.argv[1:]
    csv_file = None
    #step one: parse args
    dictionary_args = get_split_arguments(args)
    try:
        #step two: read molecule from file
        mol=MoleculeReader.from_file(**dictionary_args)
        dictionary_args['molecule']=mol
        #step three: print molecule printouts
        mol.print_equivalence_class_summary(dictionary_args['use_chains'])
        #step five: call the calculation
        if dictionary_args['calc_type'] == 'approx':
            if dictionary_args['print_approx']:
                class PrintApprox(Approx):
                    def log(self, *args, **kwargs):
                        print(*args)
                calc=PrintApprox(**dictionary_args)
            else:
                calc = Approx(**dictionary_args)
        elif dictionary_args['calc_type'] == 'trivial':
            calc = Trivial(**dictionary_args)
        else:
            try:
                # step four: create a callback function for the calculation
                # Outputing permutations
                csm_state_tracer_func= None
                if dictionary_args['perms_csv_name']:
                    csv_file = open(dictionary_args['perms_csv_name'], 'w')
                    perm_writer = csv.writer(csv_file, lineterminator='\n')
                    perm_writer.writerow(['Permutation', 'Direction', 'CSM'])
                    csm_state_tracer_func = lambda state: perm_writer.writerow(
                        [[p + 1 for p in state.perm],
                         state.dir,
                         state.csm, ])
                calc=Exact(**dictionary_args, callback_func=csm_state_tracer_func)
            except TimeoutError:
                print("Timed out")
                return
        result=calc.result
        #step six: print the results
        if dictionary_args['calc_local']:
            result.compute_local_csm()
        r=FileWriter(result, **dictionary_args)
        r.write()
        return result
    except Exception as e:
        if dictionary_args['json_output']:
            json_dict = {
                "Error": str(e)
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
    print("Runtime:", time, "seconds")
