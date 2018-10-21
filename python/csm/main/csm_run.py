import sys
import timeit
from csm.calculations.basic_calculations import CalculationTimeoutError
from csm.input_output.arguments import get_split_arguments
from csm.calculations import Approx, Trivial, Exact
from csm.input_output.writers import FileWriter
from csm import __version__
from csm.molecule.molecule import MoleculeReader

sys.setrecursionlimit(10000)

def run(args=[]):
    print("Protein-CSM version %s" % __version__)
    if not args:
        args = sys.argv[1:]

    #parse args
    dictionary_args = get_split_arguments(args)

    #read molecule from file
    mol=MoleculeReader.from_file(**dictionary_args)
    dictionary_args['molecule']=mol

    #print molecule printouts
    mol.print_equivalence_class_summary(dictionary_args['use_chains'])

    #call the calculation
    if dictionary_args['calc_type'] == 'approx':
        calc = Approx(**dictionary_args)
    elif dictionary_args['calc_type'] == 'trivial':
        calc = Trivial(**dictionary_args)
    else:
        csm_state_tracer_func= None
        calc=Exact(**dictionary_args, callback_func=csm_state_tracer_func)

    try:
        calc.calculate()
    except CalculationTimeoutError as e:
        print("Timed out")
        return

    result=calc.result
    #print the results
    fw=FileWriter(result, format=dictionary_args['molecule'].metadata.format, **dictionary_args)
    fw.write()
    return result

def run_no_return(args=[]):
    run(args)

if __name__ == '__main__':
    timer = timeit.Timer(lambda: run(args=sys.argv[1:]))
    time = timer.timeit(number=1)
    print("Runtime:", time, "seconds")
