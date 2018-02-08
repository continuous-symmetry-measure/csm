import csv
import json
import logging
import sys
import timeit

from csm.calculations.approx.approximators import ParallelApprox
from csm.calculations.constants import CalculationTimeoutError
from csm.calculations.data_classes import CSMResult, Operation
from csm.input_output.arguments import get_parsed_args
from csm.calculations import Approx, Trivial, Exact, DirectionChooser
from csm.input_output.readers import read_perm
from csm.input_output.writers import FileWriter, ApproxStatisticWriter
from csm import __version__
from csm.molecule.molecule import MoleculeReader, Molecule
from csm.input_output.formatters import csm_log as print
sys.setrecursionlimit(10000)

def read_molecule(dictionary_args):
    mol = MoleculeReader.from_file(**dictionary_args)
    mol.print_equivalence_class_summary(dictionary_args['use_chains'])
    return mol

def write_results(dictionary_args, result):
    # step six: print the results
    if "op_name" not in dictionary_args:
        op=Operation.placeholder(result.op_type, result.op_order)
        dictionary_args["op_name"]=op.name
    if dictionary_args['calc_local']:
        result.compute_local_csm()
    fw = FileWriter(result, **dictionary_args)
    fw.write()


def run(args=[]):
    print("CSM version %s" % __version__)
    if not args:
        args = sys.argv[1:]
    dictionary_args=get_parsed_args(args)

    command= dictionary_args["command"]

    if command=="read":
        mol=read_molecule(dictionary_args)
        sys.stdout.write(json.dumps(mol.to_dict(), indent=4))
    elif command == "write":
        raw_json = sys.stdin.read()
        result_dict=json.loads(raw_json)
        result=CSMResult.from_dict(result_dict)
        write_results(dictionary_args, result)
    else:
        #get input:
        if dictionary_args["in_file_name"]:
            dictionary_args['molecule']=read_molecule(dictionary_args)
        else:
            print("reading from stdin")
            raw_json=sys.stdin.read()
            print("raw json", raw_json)
            dictionary_args['molecule']=Molecule.from_json(raw_json)

        if command=="exact":
            #get perm if it exists:
            dictionary_args['perm'] = read_perm(**dictionary_args)
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

        if command=="approx":
            dir_choose = DirectionChooser(**dictionary_args)
            dictionary_args["direction_chooser"] = dir_choose
            if dictionary_args["parallel"]:
                calc=ParallelApprox(**dictionary_args)
            else:
                if dictionary_args['print_approx']:
                    def log(self, *args, **kwargs):
                        print(*args)
                    dictionary_args["log_func"]=log
                calc = Approx(**dictionary_args)

        if command=="trivial":
            calc = Trivial(**dictionary_args)

        #run the calculation
        try:
            calc.calculate()
        except CalculationTimeoutError as e:
            print("Timed out")
            return
        result=calc.result

        #statistics for exact:
        try:
            if dictionary_args["print_branches"]:
                calc.statistics.write_to_screen()
        except KeyError:
            pass

        #statistics for approx
        try:
            if dictionary_args["stat_file_name"] is not None:
                sw=ApproxStatisticWriter(calc.statistics, dictionary_args["stat_file_name"], dictionary_args["polar"])
                sw.write()
        except KeyError:
            pass

        #do output:
        if dictionary_args["out_file_name"]:
            write_results(dictionary_args, result)
        else:
            sys.stdout.write(json.dumps(result.to_dict(), indent=4))



def run_no_return(args=[]):
    run(args)


if __name__ == '__main__':
    #timer = timeit.Timer(lambda: run_many_times())
    timer = timeit.Timer(lambda: run(args=sys.argv[1:]))
    time = timer.timeit(number=1)
    print("Runtime:", time, "seconds")
