import csv
import json
import logging
import sys
import timeit
import os
from csm.calculations.approx.dirs import get_direction_chooser
from csm.calculations.constants import CalculationTimeoutError
from csm.calculations.data_classes import CSMResult, Operation
from csm.input_output.arguments import get_parsed_args
from csm.calculations import Approx, Trivial, Exact, ParallelApprox, DirectionChooser
from csm.input_output.readers import read_perm, read_from_sys_std_in
from csm.input_output.writers import OldFormatFileWriter, ScriptWriter
from csm import __version__
from csm.molecule.molecule import MoleculeReader, Molecule
from csm.input_output.formatters import csm_log as print
from csm.main.normcsm import norm_calc
sys.setrecursionlimit(10000)

def read_molecules(**kwargs):
    input_name=kwargs["in_file_name"]
    if os.path.isdir(input_name):
        x = kwargs.pop('in_file_name')
        mols=[]
        for directory, subdirectories, files in os.walk(input_name):
            for file_name in files:
                if file_name=="sym.txt":
                    continue
                mol_file= os.path.join(input_name, file_name)
                try:
                    mol = MoleculeReader.from_file(mol_file, **kwargs)
                    mols.append(mol)
                except Exception as e:
                    print("failed to create a molecule from", file_name, e)
    elif not os.path.isfile(input_name):
        raise ValueError("invalid file/folder name for molecule")
    else: #file
        mols=MoleculeReader.multiple_from_file(**kwargs)
    sys.stderr.flush()

    if kwargs['select_mols']:
        try:
            mols=[mols[i] for i in kwargs['select_mols']]
        except IndexError:
            raise IndexError("You have selected more molecules than you have input")
    for mol in mols:
        mol.print_equivalence_class_summary(True)
    return mols

def write_results(results_arr, **kwargs):
    if kwargs['simple']:
        for mol_index, mol_result in enumerate(results_arr):
            for lin_index, line_result in enumerate(mol_result):
                print("mol", mol_index, "cmd", lin_index, " CSM: ", line_result.csm)
        return

    if kwargs['out_file_name']:
        if kwargs['legacy']:
            if len(results_arr)>1 or len(results_arr[0])>1:
                raise ValueError("Legacy result writing only works for a single molecule and single command")
            result=results_arr[0][0]
            writer=OldFormatFileWriter(result, **kwargs)
            writer.write()
            return


        if not os.path.isdir(kwargs['out_file_name']):
            if len(results_arr) == 1 and len(results_arr[0]) == 1:
                print("You are running a single file and command. Did you want to print to the old format, --legacy?")
            kwargs['out_file_name']=os.path.dirname(kwargs['out_file_name'])

        if 'out_format' in kwargs and kwargs['out_format']:
            format=kwargs['out_format']
        elif 'in_format' in kwargs and kwargs['in_format']:
            format=kwargs['in_format']
        else:
            format=results_arr[0][0].molecule._format
        writer = ScriptWriter(results_arr, format, **kwargs)
        writer.write()
        return

    #default option
    sys.stdout.write(json.dumps([[result.to_dict() for result in mol_results_arr] for mol_results_arr in results_arr], indent=4))





def do_calculation(dictionary_args):
    command = dictionary_args["command"]
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
        dir_chooser = get_direction_chooser(**dictionary_args)
        dictionary_args["direction_chooser"] = dir_chooser
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
    return calc.result


def run(args=[]):
    print("CSM version %s" % __version__)
    #get command
    if not args:
        args = sys.argv[1:]
    dictionary_args=get_parsed_args(args)
    command= dictionary_args["command"]

    #call command funcs:
    if command=="read":
        mols= read_molecules(**dictionary_args)
        sys.stdout.write(json.dumps([mol.to_dict() for mol in mols], indent=4))
        return mols

    elif command=="write":
        raw_json = read_from_sys_std_in()
        less_raw_json = json.loads(raw_json)
        results=[[CSMResult.from_dict(result_dict) for result_dict in mol_arr] for mol_arr in less_raw_json]
        write_results(results, **dictionary_args)

    else:
        if dictionary_args["in_file_name"]:
            molecules = read_molecules(**dictionary_args)
        else:
            raw_json = read_from_sys_std_in()
            less_raw_json=json.loads(raw_json)
            molecules=[Molecule.from_dict(json_dict) for json_dict in less_raw_json]

        if command == "command":
            args_array = []
            command_file = dictionary_args["command_file"]
            with open(command_file, 'r') as file:
                for line in file:
                    try:
                        cmd_arg = get_parsed_args(line.split())
                        args_array.append(cmd_arg)
                    except: #want to be able to run even if some lines are invalid
                        print("failed to read args from line", line)
        else:
            args_array=[dictionary_args]

        total_results=[]

        for molecule in molecules:
            mol_results=[]
            for line_args in args_array:
                line_args['molecule'] = molecule
                result = do_calculation(line_args)
                mol_results.append(result)
                try:
                    if len(dictionary_args['normalizations']) > 0:
                        norm_calc(result, dictionary_args['normalizations'])
                except KeyError:
                    pass
            total_results.append(mol_results)
        write_results(total_results, **dictionary_args)
        return total_results

def run_no_return(args=[]):
    run(args)


if __name__ == '__main__':
    timer = timeit.Timer(lambda: run(args=sys.argv[1:]))
    time = timer.timeit(number=1)
    print("Runtime:", time, "seconds")
