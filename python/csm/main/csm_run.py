import json
import logging
import sys
import timeit

from shutil import copyfile

import os

from csm.calculations.constants import CalculationTimeoutError
from csm.calculations.data_classes import CSMResult, Operation
from csm.input_output.arguments import get_parsed_args, old_cmd_converter, check_modifies_molecule
from csm import __version__
from csm.input_output.readers import read_from_sys_std_in
from csm.main.calculate import do_calculation
from csm.main.read import read_molecules, read_mols_from_std_in, read
from csm.main.write import write_results, write

from csm.input_output.formatters import csm_log as print
from csm.main.normcsm import norm_calc
from csm.molecule.molecule import MoleculeReader

sys.setrecursionlimit(10000)


def sym_run():
    #by default, sym_run reads any molecules in the folder it's being run from

    #by default, it reads the command.txt in that folder

    pass


def get_command_args(command_file, old_command=True):
    args_array=[]
    with open(command_file, 'r') as file:
        for line in file:
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
            args_array.append((args_dict, modifies_molecule))
    return args_array

def do_commands(molecules, **dictionary_args):
    if not os.path.isdir(dictionary_args["out_file_name"]):
        os.mkdir(dictionary_args["out_file_name"])
    copyfile(dictionary_args["command_file"], os.path.join(dictionary_args["out_file_name"], "command.txt"))


    args_array=get_command_args(dictionary_args["command_file"], dictionary_args["old_command"])
    total_results=[[] for mol in molecules]
    for args_dict, modifies_molecule in args_array:
            for mol_index, molecule in enumerate(molecules):
                args_dict["molecule"]=molecule
                if modifies_molecule:
                    new_molecule=MoleculeReader.redo_molecule(molecule, **args_dict)
                    args_dict["molecule"]=new_molecule
                result = do_calculation(**args_dict)
                total_results[mol_index].append(result)

    write_results(total_results, **dictionary_args)
    return total_results

def csm_run(args=[]):
    print("CSM version %s" % __version__)
    #get command
    if not args:
        args = sys.argv[1:]
    dictionary_args=get_parsed_args(args)
    command= dictionary_args["command"]

    #call command funcs:
    if command=="read":
        read(**dictionary_args)

    elif command=="write":
        write(**dictionary_args)

    else:
        if dictionary_args["in_file_name"]:
            molecules = read_molecules(**dictionary_args)
        else:
            molecules = read_mols_from_std_in()

        if command=="command":
            return do_commands(molecules, **dictionary_args)

        total_results=[]
        for molecule in molecules:
            dictionary_args["molecule"]=molecule
            try:
                result = do_calculation(**dictionary_args)
                total_results.append(result)
                try:
                    if len(dictionary_args['normalizations']) > 0:
                        norm_calc(result, dictionary_args['normalizations'])
                except KeyError:
                    pass
            except CalculationTimeoutError as e:
                print("Calculation timed out")

        write_results(total_results, **dictionary_args)
        return total_results

def run_no_return(args=[]):
    csm_run(args)


if __name__ == '__main__':
    timer = timeit.Timer(lambda: csm_run(args=sys.argv[1:]))
    time = timer.timeit(number=1)
    print("Runtime:", time, "seconds")
