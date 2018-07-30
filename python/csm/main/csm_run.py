import sys
import timeit

from shutil import copyfile

import os

from csm.calculations.constants import CalculationTimeoutError
from csm.calculations.data_classes import FailedResult
from csm.input_output.arguments import get_parsed_args, old_cmd_converter, check_modifies_molecule
from csm import __version__
from csm.main.calculate import single_calculation
from csm.input_output.read import read_molecules, read_mols_from_std_in, read
from csm.input_output.write import write_results, write

from csm.input_output.formatters import csm_log as print
from csm.main.normcsm import norm_calc
from csm.molecule.molecule import MoleculeReader

sys.setrecursionlimit(10000)

def get_command_args(command_file, old_command=True):
    args_array=[]
    with open(command_file, 'r') as file:
        for line in file:
            if line[0]=="#":
                continue
            modifies_molecule=check_modifies_molecule(line)
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
    return args_array

def do_commands(molecules, **dictionary_args):
    if not os.path.isdir(dictionary_args["out_file_name"]):
        os.makedirs(dictionary_args["out_file_name"], exist_ok=True)
    copyfile(dictionary_args["command_file"], os.path.join(dictionary_args["out_file_name"], "command.txt"))


    args_array=get_command_args(dictionary_args["command_file"], dictionary_args["old_command"])
    total_results=[[] for mol in molecules]
    for line, args_dict, modifies_molecule in args_array:
        print("\nexecuting command:", line[:-1])
        try:
            selections=args_dict['select_mols']
            actual_mols = [molecules[i] for i in selections]
            assert len(actual_mols)>0
        except IndexError:
            raise IndexError("You have selected more molecules than you have input")
        except AssertionError:
            actual_mols = molecules

        for mol_index, molecule in enumerate(actual_mols):
            molecule.print_equivalence_class_summary(True)
            args_dict["molecule"]=molecule
            if modifies_molecule:
                new_molecule=MoleculeReader.redo_molecule(molecule, **args_dict)
                new_molecule.metadata.index=mol_index
                args_dict["molecule"]=new_molecule
            try:
                result= single_calculation(args_dict["molecule"], args_dict)
                total_results[mol_index].append(result)
            except Exception as e:#CalculationTimeoutError:
                print(str(e))
                total_results[mol_index].append(FailedResult(str(e), **dictionary_args))

    return total_results


def csm_run(args=[]):
    #get command
    if not args:
        args = sys.argv[1:]
    print(" ".join(args))
    dictionary_args=get_parsed_args(args)
    if dictionary_args["pipe"]:
        from csm.input_output import formatters
        formatters.csm_out_pipe=sys.stderr

    print("CSM version %s" % __version__)
    command= dictionary_args["command"]

    #call command funcs:
    if command=="read":
        return read(**dictionary_args)

    elif command=="write":
        return write(**dictionary_args)


    if dictionary_args["in_file_name"]:
        molecules = read_molecules(**dictionary_args)
    elif dictionary_args["pipe"]:
        molecules = read_mols_from_std_in()
    else:
        raise ValueError("No input for molecules specified")

    if True:
        if command=="command":
            total_results= do_commands(molecules, **dictionary_args)
        else:
            total_results=[]
            for molecule in molecules:
                try:
                    result= single_calculation(molecule, dictionary_args)
                    total_results.append([result])
                except Exception as e:#CalculationTimeoutError:
                    print(str(e))
                    total_results.append([FailedResult(str(e), **dictionary_args)])

        write_results(total_results, **dictionary_args)
        return total_results



def run_no_return(args=[]):
    csm_run(args)


if __name__ == '__main__':
    timer = timeit.Timer(lambda: csm_run(args=sys.argv[1:]))
    time = timer.timeit(number=1)
    print("Runtime: "+ str(time)+ " seconds")
