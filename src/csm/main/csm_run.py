import csv
import json
import multiprocessing
import os
import sys
import timeit
import numpy as np

sys.path.append('src')

from csm.main.openbabel_fix import prepare_openbabel
prepare_openbabel()  # Make sure the OpenBabel DLL can be found. See file for more information

from csm import __version__
from csm.input_output import formatters
from csm.calculations import Approx, Trivial, Exact, ParallelApprox
from csm.calculations.approx.dirs import get_direction_chooser
from csm.calculations.data_classes import FailedResult, CSMResult
from csm.input_output.arguments import get_parsed_args, old_cmd_converter, check_modifies_molecule
from csm.input_output.formatters import csm_log as print
from csm.input_output.formatters import silent_print
from csm.input_output.readers import read_molecules, read_mols_from_std_in, read
from csm.input_output.readers import read_perm, read_from_sys_std_in
from csm.input_output.writers import SimpleContextWriter, ScriptContextWriter, PipeContextWriter, LegacyContextWriter, \
    get_line_header, MoleculeWriter
from csm.molecule.molecule import Molecule
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
                [state.op_type + str(state.op_order),
                 [p + 1 for p in state.perm],
                 state.dir,
                 state.csm, ])
        calc = Exact(**dictionary_args, callback_func=csm_state_tracer_func)

    elif calc_type == "approx":
        dictionary_args['chain_perms'] = read_perm(**dictionary_args)
        dir_chooser = get_direction_chooser(**dictionary_args)
        dictionary_args["direction_chooser"] = dir_chooser
        if parallel_dirs:
            parallel_obmol = dictionary_args["molecule"]._obmol
            dictionary_args["molecule"]._obmol = None
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
    if parallel_dirs:
        # manage pickling
        dictionary_args["molecule"]._obmol = parallel_obmol
        calc.result.molecule._obmol = parallel_obmol

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
                args_dict = {**dictionary_args, **in_args}
                args_dict["line_command"] = line
            except:  # want to be able to run even if some lines are invalid
                print("failed to read args from line", line)
                continue
            if args_dict.get('select_mols', False):
                raise ValueError(
                    "Not allowed to use argument --select-mols within command file. Please use it in the main program command (eg `comfile cmd.txt --select-mols 1-3`)")
            args_array.append((line, args_dict, modifies_molecule))
            operation_array.append(args_dict["operation"])
    return args_array, operation_array


def csm_run(args=[]):
    formatters.csm_out_pipe = sys.stdout
    # get command
    if not args:
        args = sys.argv[1:]
    if 'read' not in args:
        print("CSM version %s" % __version__)
        print(" ".join(args))

    dictionary_args = get_parsed_args(args)
    dictionary_args["argument_string"] = " ".join(args) + "\n"
    if "global_timeout" in dictionary_args:
        from csm.calculations.constants import set_global_timeout
        set_global_timeout(dictionary_args["global_timeout"])
    if dictionary_args["pipe"]:
        formatters.csm_out_pipe = sys.stderr

    command = dictionary_args["command"]

    # call command func that aren't calculate:
    if command == "read":
        formatters.csm_out_pipe = sys.stderr
        try:
            return read(**dictionary_args)
        except Exception as ex:
            formatters.csm_out_pipe = sys.stdout
            print("error in read", str(ex))
            return {"error in read": str(ex)}

    elif command == "write":
        return write(**dictionary_args)

    else:
        try:
            return calc(dictionary_args)
        except Exception as err:
            print(err)
            exit(2)


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
    if "error in read" in raw_json:
        print("Error in the read step:", raw_json.rsplit('error in read', 1)[-1], file=sys.stdout)
        return
    less_raw_json = json.loads(raw_json)
    
    mols: list[Molecule] = [Molecule.from_dict(mol_dict) for mol_dict in less_raw_json]
    if not mols:
        raise ValueError("write- Not found molecules")

    out_format = mols[0].metadata.format
    out_filename = dictionary_args.get('out_file_name', 'result-molecule')
    if not out_filename.endswith(out_format):
        out_filename = out_filename + f'.{out_format}'
    mol_writer = MoleculeWriter(None, out_format=out_format)
    model_number = 0
    with open(out_filename, 'w') as file:
        for molecule in mols:
            mol_writer.write(file, np.array(molecule.Q), consecutive=True, model_number=model_number, obmols=molecule.obmol)
            model_number += 1
        if out_format.lower() == 'pdb':
            file.write("\nEND\n")
    print(f"The result saved on the file: {out_filename}")




def calc(dictionary_args):
    # get commands:
    if dictionary_args["command"] == "comfile":
        args_array, operation_array = get_command_args(**dictionary_args)
    else:
        args_array = [(None, dictionary_args, False)]
        operation_array = [dictionary_args["operation"]]

    # get molecules
    if dictionary_args["in_file_name"]:
        dictionary_args['comfile_first_read'] = dictionary_args["command"] == "comfile"
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
        try:
            mol_args = []
            for line, args_dict, modifies_molecule in args_array:
                args_dict["molecule"] = molecule
                args_dict["line"] = line
                # handle modifying molecules:
                if modifies_molecule:
                    new_molecule = MoleculeReader.redo_molecule(molecule, **args_dict)
                    new_molecule.metadata.index = mol_index
                    args_dict["molecule"] = new_molecule
                mol_args.append(dict(args_dict))
            total_args.append(mol_args)
        except Exception as ex:
            print(f"Error with the molecule {molecule.metadata.filepath}:")
            print(ex)
            pass # continue with the rest of the molecules

    # run the calculation, in parallel
    if dictionary_args["parallel"]:
        flattened_args = [item for sublist in total_args for item in sublist]
        # manage pickling
        flattened_args_obmols = [dic_arg["molecule"].obmol for dic_arg in flattened_args]
        for dic_arg in flattened_args:
            dic_arg["molecule"]._obmol = None

        num_ops = len(operation_array)
        batch_mols = 50  # int(len(molecules)/10)
        batch_size = num_ops * batch_mols  # it needs to be divisible by length of operation array
        total_results = []
        pool_size = dictionary_args["pool_size"]
        print("Parallelizing {} calculations across {} processes with batch size {}".format(len(flattened_args),
                                                                                            pool_size, batch_size))
        try:
            pool = multiprocessing.Pool(processes=pool_size)
            with context_writer(operation_array, **dictionary_args) as rw:
                for i in range(0, len(flattened_args), batch_size):
                    end_index = min(i + batch_size, len(flattened_args))
                    args_array = flattened_args[i:end_index]
                    now = datetime.now()
                    # print("calculating partial results for chunk{}-{} - {}".format(i, end_index, now.strftime("%d/%m/%Y %H:%M:%S")))
                    partial_results = pool.map(single_calculation, args_array)
                    obmol_array = flattened_args_obmols[i:end_index]

                    for index in range(len(obmol_array)):
                        partial_results[index].molecule._obmol = obmol_array[index]

                    m_range = int((end_index - i) / num_ops)
                    unflattened_partial_results = [
                        [partial_results[m_index * num_ops + command_index] for command_index in range(num_ops)] for
                        m_index in
                        range(m_range)
                    ]
                    # now=datetime.now()
                    # print("outputting partial results for chunk{}-{} - {}".format(i, end_index, now.strftime("%d/%m/%Y %H:%M:%S")))
                    for mol_results in unflattened_partial_results:
                        rw.write(mol_results)
                    total_results = total_results + partial_results
        except Exception as e:
            print(e)
        finally:
            pool.close()
            pool.join()
        return total_results  # maybe should unflatten first?

    if len(molecules) > 10:
        from csm.input_output.formatters import output_strings
        output_strings.silent = True
        print(len(molecules),
              " molecules in folder. Molecule and result summaries can be found in extra.txt and will not be printed to screen")

    # run the calculation, in serial
    all_results = []
    with context_writer(operation_array,
                        **dictionary_args) as rw:  # this is the line of code where the results folder is created
        for mol_index, mol_args in enumerate(total_args):
            mol_results = []
            for line_index, args_dict in enumerate(mol_args):

                # create perms.csv if relevant
                if args_dict.get('output_perms', False):
                    args_dict["perms_csv_name"] = rw.create_perms_csv(args_dict, line_index)

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


def run_no_return(args=[]):
    csm_run(args)


if __name__ == '__main__':
    timer = timeit.Timer(lambda: csm_run(args=sys.argv[1:]))
    time = timer.timeit(number=1)
    if 'read' not in sys.argv[1:]:
        print("-----\nRuntime: " + str(time) + " seconds")
