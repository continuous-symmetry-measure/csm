import csv
import json
import multiprocessing
import os
import sys
import timeit

from csm import __version__
from csm.calculations import Approx, Trivial, Exact, ParallelApprox
from csm.calculations.approx.dirs import get_direction_chooser
from csm.calculations.data_classes import FailedResult, CSMResult
from csm.input_output.arguments import get_parsed_args, old_cmd_converter, check_modifies_molecule
from csm.input_output.formatters import csm_log as print
from csm.input_output.formatters import silent_print
from csm.input_output.readers import read_molecules, read_mols_from_std_in, read
from csm.input_output.readers import read_perm, read_from_sys_std_in
from csm.input_output.writers import SimpleContextWriter, ScriptContextWriter, PipeContextWriter, LegacyContextWriter, \
    get_line_header
from csm.main.normcsm import norm_calc
from csm.molecule.molecule import MoleculeReader
from datetime import datetime

def do_calculation(command, perms_csv_name=None, parallel_dirs=False, print_approx=False, **dictionary_args):
    calc_type = command
    if calc_type == "exact":
        # get perm if it exists:
        dictionary_args['perm'] = read_perm(**dictionary_args)
        csm_state_tracer_func = None
        if perms_csv_name:
            csv_file = open(perms_csv_name, 'a')
            perm_writer = csv.writer(csv_file, lineterminator='\n')
            csm_state_tracer_func = lambda state: perm_writer.writerow(
                [state.op_type+str(state.op_order),
                [p + 1 for p in state.perm],
                 state.dir,
                 state.csm, ])
        calc = Exact(**dictionary_args, callback_func=csm_state_tracer_func)

    elif calc_type == "approx":
        dictionary_args['chain_perms'] = read_perm(**dictionary_args)
        dir_chooser = get_direction_chooser(**dictionary_args)
        dictionary_args["direction_chooser"] = dir_chooser
        if parallel_dirs:
            calc = ParallelApprox(**dictionary_args)
        else:
            if print_approx:
                def log(self, *args, **kwargs):
                    print(*args)

                dictionary_args["log_func"] = log
            calc = Approx(**dictionary_args)

    elif calc_type == "trivial":
        dictionary_args['chain_perms'] = read_perm(**dictionary_args)
        calc = Trivial(**dictionary_args)

    # run the calculation
    calc.calculate(**dictionary_args)
    return calc.result


def single_calculation(dictionary_args):
    result = do_calculation(**dictionary_args)
    try:
        if len(dictionary_args['normalizations']) > 0:
            norm_calc(result, dictionary_args['normalizations'])
    except KeyError:
        pass
    return result


def get_command_args(command_file, old_command=True, **dictionary_args):
    args_array = []
    operation_array = []
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
                in_args = get_parsed_args(fixed_args)
                args_dict={**dictionary_args, **in_args}
                args_dict["line_command"]=line
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
    dictionary_args["argument_string"] = " ".join(args) + "\n"
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
    elif dictionary_args["legacy_output"]:
        context_writer = LegacyContextWriter
    else:
        context_writer = ScriptContextWriter
    return context_writer


def write(**dictionary_args):
    raw_json = read_from_sys_std_in()
    less_raw_json = json.loads(raw_json)
    results = [[CSMResult.from_dict(result_dict) for result_dict in mol_arr] for mol_arr in less_raw_json]
    context_writer = get_context_writer(dictionary_args)
    writer = context_writer(results, context_writer=context_writer, **dictionary_args)
    writer.write()


def calc(dictionary_args):
    # get commands:
    if dictionary_args["command"] == "comfile":
        args_array, operation_array = get_command_args(**dictionary_args)
    else:
        args_array = [(None, dictionary_args, False)]
        operation_array = [dictionary_args["operation"]]


    # get molecules
    if dictionary_args["in_file_name"]:
        molecules = read_molecules(**dictionary_args)
    elif dictionary_args["pipe"]:
        molecules = read_mols_from_std_in()
    else:
        raise ValueError("No input for molecules specified")

    max_len_file_name = 0
    for mol in molecules:
        len_name = len(mol.metadata.filename)
        if len_name > max_len_file_name:
            max_len_file_name = len_name
    dictionary_args["max_len_file_name"] = max_len_file_name

    context_writer = get_context_writer(dictionary_args)
    if not dictionary_args["out_format"]:
        dictionary_args["out_format"] = dictionary_args["in_format"]
    if not dictionary_args["out_format"]:
        dictionary_args["out_format"] = molecules[0].metadata.format

    # process arguments into flat and unflat arrays
    total_args = []
    for mol_index, molecule in enumerate(molecules):
        mol_args = []
        for line, args_dict, modifies_molecule in args_array:
            args_dict["molecule"] = molecule
            args_dict["line"] = line
            # handle select molecules:
            if args_dict.get('select_mols', False):
                raise ValueError("Not allowed to use argument --select-mols within command file. Please use it in the main program command (eg comfile cmd.txt --select-mols)")


            # handle modifying molecules:
            if modifies_molecule:
                new_molecule = MoleculeReader.redo_molecule(molecule, **args_dict)
                new_molecule.metadata.index = mol_index
                args_dict["molecule"] = new_molecule
            mol_args.append(dict(args_dict))
        total_args.append(mol_args)

    # run the calculation, in parallel
    if dictionary_args["parallel"]:
        flattened_args = [item for sublist in total_args for item in sublist]
        num_ops = len(operation_array)
        batch_mols= 50  # int(len(molecules)/10)
        batch_size = num_ops * batch_mols #it needs to be divisible by length of operation array
        total_results = []
        pool_size=dictionary_args["pool_size"]
        print("Parallelizing {} calculations across {} processes with batch size {}".format(len(flattened_args), pool_size, batch_size))
        try:
            pool = multiprocessing.Pool(processes=pool_size)
            with context_writer(operation_array, **dictionary_args) as rw:
                for i in range(0, len(flattened_args), batch_size):
                    end_index=min(i+batch_size, len(flattened_args))
                    args_array=flattened_args[i:end_index]
                    now=datetime.now()
                    #print("calculating partial results for chunk{}-{} - {}".format(i, end_index, now.strftime("%d/%m/%Y %H:%M:%S")))
                    partial_results = pool.map(single_calculation, args_array)



                    m_range=int((end_index-i)/num_ops)
                    unflattened_partial_results=[
                    [partial_results[m_index * num_ops + command_index] for command_index in range(num_ops)] for m_index in
                    range(m_range)
                ]
                    #now=datetime.now()
                    #print("outputting partial results for chunk{}-{} - {}".format(i, end_index, now.strftime("%d/%m/%Y %H:%M:%S")))
                    for mol_results in unflattened_partial_results:
                        rw.write(mol_results)
                    total_results=total_results+partial_results
        except Exception as e:
            print(e)
        finally:
            pool.close()
            pool.join()
        return total_results #maybe should unflatten first?

    if len(molecules) > 10:
        from csm.input_output.formatters import output_strings
        output_strings.silent = True
        print(len(molecules),
              " molecules in folder. Molecule and result summaries can be found in extra.txt and will not be printed to screen")

    # run the calculation, in serial
    all_results = []
    with context_writer(operation_array, **dictionary_args) as rw:  #this is the line of code where the results folder is created
        for mol_index, mol_args in enumerate(total_args):
            mol_results = []
            for line_index, args_dict in enumerate(mol_args):
                
                # create perms.csv if relevant
                args_dict["perms_csv_name"]=create_perms_csv(args_dict, line_index, rw)

                # print stuff
                molecule = args_dict["molecule"]

                if args_dict["line"]:
                    silent_print("-----")
                    silent_print("**executing command:", args_dict["line"].rstrip(), "**")
                silent_print("Molecule:", molecule.metadata.appellation(no_leading_zeros=True))
                molecule.print_equivalence_class_summary(True)
                # run the calculation
                try:
                    result = single_calculation(args_dict)
                    result.print_summary(dictionary_args["legacy_output"])
                except Exception as e:
                    print(e)
                    result = FailedResult(str(e), **args_dict)
                mol_results.append(result)
            # write the results for the molecule
            rw.write(mol_results)
            all_results.append(mol_results)
    return all_results

def create_perms_csv(args_dict, line_index, rw):
    perms_csv_name=None
    output_perms = args_dict.get("output_perms")
    if output_perms:
        filename = args_dict["molecule"].metadata.appellation(no_file_format=True) + "_" + get_line_header(line_index,
                                                                                        args_dict["operation"]) + ".csv"
        perms_csv_name = os.path.join(rw.folder, 'exact', filename)
        csv_file = open(perms_csv_name, 'w')
        perm_writer = csv.writer(csv_file, lineterminator='\n')
        perm_writer.writerow(['op', 'Permutation', 'Direction', 'CSM'])
    return perms_csv_name

def run_no_return(args=[]):
    csm_run(args)


if __name__ == '__main__':
    timer = timeit.Timer(lambda: csm_run(args=sys.argv[1:]))
    time = timer.timeit(number=1)
    print("-----\nRuntime: " + str(time) + " seconds")
