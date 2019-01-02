import json
import sys
import timeit

import os
from shutil import copyfile

from csm import __version__
from csm.calculations.data_classes import FailedResult, CSMResult
from csm.input_output.arguments import get_parsed_args, old_cmd_converter, check_modifies_molecule
from csm.input_output.formatters import csm_log as print
from csm.input_output.read import read_molecules, read_mols_from_std_in, read
from csm.input_output.writers import SimpleContextWriter, ScriptContextWriter, PipeContextWriter, LegacyContextWriter, ConsolidatedScriptWriter
from csm.molecule.molecule import MoleculeReader
import csv

from csm.calculations import Approx, Trivial, Exact, ParallelApprox
from csm.calculations.approx.dirs import get_direction_chooser
from csm.input_output.formatters import csm_log as print
from csm.input_output.readers import read_perm, read_from_sys_std_in
from csm.main.normcsm import norm_calc

def do_calculation(command, perms_csv_name=None, parallel=False, print_approx=False, **dictionary_args):
    calc_type = command
    if calc_type == "exact":
        # get perm if it exists:
        dictionary_args['perm'] = read_perm(**dictionary_args)
        csm_state_tracer_func = None
        if perms_csv_name:
            csv_file = open(perms_csv_name, 'w')
            perm_writer = csv.writer(csv_file, lineterminator='\n')
            perm_writer.writerow(['Permutation', 'Direction', 'CSM'])
            csm_state_tracer_func = lambda state: perm_writer.writerow(
                [[p + 1 for p in state.perm],
                 state.dir,
                 state.csm, ])
        calc = Exact(**dictionary_args, callback_func=csm_state_tracer_func)

    if calc_type == "approx":
        dir_chooser = get_direction_chooser(**dictionary_args)
        dictionary_args["direction_chooser"] = dir_chooser
        if parallel:
            calc = ParallelApprox(**dictionary_args)
        else:
            if print_approx:
                def log(self, *args, **kwargs):
                    print(*args)

                dictionary_args["log_func"] = log
            calc = Approx(**dictionary_args)

    if calc_type == "trivial":
        calc = Trivial(**dictionary_args)

    # run the calculation
    calc.calculate(**dictionary_args)
    return calc.result


def single_calculation(dictionary_args):
    molecule=dictionary_args["molecule"]
    print("Molecule:", molecule.metadata.appellation(no_leading_zeros=True))
    molecule.print_equivalence_class_summary(True)
    dictionary_args["molecule"] = molecule
    result = do_calculation(**dictionary_args)
    result.print_summary()
    try:
        if len(dictionary_args['normalizations']) > 0:
            norm_calc(result, dictionary_args['normalizations'])
    except KeyError:
        pass
    print("-----")
    return result

def get_command_args(command_file, old_command=True):
    args_array = []
    operation_array=[]
    with open(command_file, 'r') as file:
        for line in file:
            if line[0] == "#":
                continue
            modifies_molecule = check_modifies_molecule(line)
            if old_command:
                fixed_args = old_cmd_converter(line)
            else:
                fixed_args = line.split()
            try:
                args_dict = get_parsed_args(fixed_args)
            except:  # want to be able to run even if some lines are invalid
                print("failed to read args from line", line)
                continue
            args_array.append((line, args_dict, modifies_molecule))
            operation_array.append(args_dict["operation"])
    return args_array, operation_array

def csm_run(args=[]):
  # get command
    if not args:
        args = sys.argv[1:]
    print("CSM version %s" % __version__)
    print(" ".join(args))

    dictionary_args = get_parsed_args(args)
    dictionary_args["argument_string"]= " ".join(args) + "\n"
    if "global_timeout" in dictionary_args:
        from csm.calculations.constants import set_global_timeout
        set_global_timeout(dictionary_args["global_timeout"])
    if dictionary_args["pipe"]:
        from csm.input_output import formatters
        formatters.csm_out_pipe = sys.stderr

    command = dictionary_args["command"]

    # call command funcs that aren't calculate:
    if command == "read":
        return read(**dictionary_args)

    elif command == "write":
        return write(**dictionary_args)

    else:
        return calc(dictionary_args)


def get_context_writer(dictionary_args):
    if dictionary_args["simple"]:
        context_writer = SimpleContextWriter
    elif dictionary_args["pipe"]:
        context_writer = PipeContextWriter
    elif dictionary_args["legacy"]:
        context_writer = LegacyContextWriter
    else:
        context_writer = ScriptContextWriter
    return context_writer

def write(**dictionary_args):
    raw_json = read_from_sys_std_in()
    less_raw_json = json.loads(raw_json)
    results = [[CSMResult.from_dict(result_dict) for result_dict in mol_arr] for mol_arr in less_raw_json]
    context_writer=get_context_writer(dictionary_args)
    writer = ConsolidatedScriptWriter(results, context_writer=context_writer, **dictionary_args)
    writer.write()

def calc(dictionary_args):
    # get molecules
    if dictionary_args["in_file_name"]:
        molecules = read_molecules(**dictionary_args)
    elif dictionary_args["pipe"]:
        molecules = read_mols_from_std_in()
    else:
        raise ValueError("No input for molecules specified")

    print("----------")

    #get commands:
    if dictionary_args["command"] == "comfile":
        args_array, operation_array = get_command_args(dictionary_args["command_file"], dictionary_args["old_command"])
    else:
        args_array = [(None, dictionary_args, False)]
        operation_array=[dictionary_args["operation"]]

    context_writer = get_context_writer(dictionary_args)
    if not dictionary_args["out_format"]:
        dictionary_args["out_format"]=dictionary_args["in_format"]
    if not dictionary_args["out_format"]:
        dictionary_args["out_format"] = molecules[0].metadata.format

    with context_writer(operation_array, **dictionary_args) as rw:
        for mol_index, molecule in enumerate(molecules):
            mol_results=[]
            for line, args_dict, modifies_molecule in args_array:
                args_dict["molecule"]=molecule
                if line:
                    print("\n**executing command:", line.rstrip(), "**")

                #handle select molecules:
                selections = args_dict['select_mols']
                if len(selections)>0 and mol_index not in selections:
                    mol_results.append(FailedResult("molecule not selected", molecule))
                    continue
                #handle modifying molecules:
                if modifies_molecule:
                    new_molecule = MoleculeReader.redo_molecule(molecule, **args_dict)
                    new_molecule.metadata.index = mol_index
                    args_dict["molecule"] = new_molecule

                #run the calculation
                try:
                    result = single_calculation(args_dict)
                except Exception as e:
                    result=FailedResult(str(e), args_dict["molecule"])
                mol_results.append(result)
            #write the results for the molecule
            rw.write(mol_results)





def run_no_return(args=[]):
    csm_run(args)


if __name__ == '__main__':
    timer = timeit.Timer(lambda: csm_run(args=sys.argv[1:]))
    time = timer.timeit(number=1)
    print("Runtime: " + str(time) + " seconds")