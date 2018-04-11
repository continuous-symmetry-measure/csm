import csv
import json
import logging
import sys
import timeit

from csm.calculations.approx.dirs import get_direction_chooser
from csm.calculations.constants import CalculationTimeoutError
from csm.calculations.data_classes import CSMResult, Operation
from csm.input_output.arguments import get_parsed_args
from csm.calculations import Approx, Trivial, Exact, ParallelApprox, DirectionChooser
from csm.input_output.readers import read_perm, read_from_sys_std_in
from csm.input_output.writers import OldFormatFileWriter, ApproxStatisticWriter
from csm import __version__
from csm.molecule.molecule import MoleculeReader, Molecule
from csm.input_output.formatters import csm_log as print
from csm.main.normcsm import norm_calc
sys.setrecursionlimit(10000)

def read_molecule(dictionary_args):
    mol = MoleculeReader.from_file(**dictionary_args)
    mol.print_equivalence_class_summary(dictionary_args['use_chains'])
    return mol

def write_results(dictionary_args, result):
    # step six: print the results
    if dictionary_args['calc_local']:
        result.compute_local_csm()
    fw = OldFormatFileWriter(result, **dictionary_args)
    fw.write()


def run_calculation(dictionary_args):
    #get input:
    if dictionary_args["in_file_name"]:
        dictionary_args['molecule']=read_molecule(dictionary_args)
    else:
        raw_json=read_from_sys_std_in()
        dictionary_args['molecule']=Molecule.from_json(raw_json)
    calc, result, dictionary_args=_run_calculation(dictionary_args)
    return _do_output(calc, result, dictionary_args)

def _run_calculation(dictionary_args):
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
    result=calc.result
    return calc, result, dictionary_args

def _do_output(calc, result, dictionary_args):
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


    if len(dictionary_args['normalizations'])>0:
        norm_calc(result, dictionary_args['normalizations'])

    return result

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
        raw_json = read_from_sys_std_in()
        result_dict=json.loads(raw_json)
        result=CSMResult.from_dict(result_dict)
        write_results(dictionary_args, result)
    else:
        return run_calculation(dictionary_args)


def folder_script_runner():
    import os
    from csm.input_output.writers import ScriptWriter
    from csm.molecule.molecule import  get_format
    in_folder=r"C:\Users\devora\Sources\temp\csm_multi_test\input"
    out_folder=r"C:\Users\devora\Sources\temp\csm_multi_test\output"

    def get_molecules(in_folder):
        molecules = []
        mol_names=[]
        for directory, subdirectories, files in os.walk(in_folder):
            for file_name in files:
                if file_name=="sym.txt":
                    continue
                mol_file= os.path.join(in_folder, file_name)
                mol = MoleculeReader.from_file(mol_file)
                mol_names.append(file_name)
                molecules.append(mol)

        return molecules, mol_names

    def get_commands(in_folder):
        commands = []
        command_file=os.path.join(in_folder, "sym.txt")
        with open(command_file, 'r') as file:
            for line in file:
                commands.append(line)
        return commands

    molecules, mol_names=get_molecules(in_folder)

    commands=get_commands(in_folder)
    dictionary_args_array=[]
    for command in commands:
        args=command.split()
        dictionary_args = get_parsed_args(args)
        dictionary_args_array.append(dictionary_args)

    results=[]
    for name, mol in zip(mol_names, molecules):
        mol_results=[]
        for args in dictionary_args_array:
            args['molecule']=mol
            calc, result, dictionary_args=_run_calculation(args)
            mol_results.append(result)
        results.append((name, mol_results))
    format=get_format(None, filename=mol_names[0])
    writer=ScriptWriter(results, format, commands, out_folder)
    writer.write()


def molecule_script_runner():
    import os
    from csm.molecule.molecule import  get_format
    from csm.input_output.writers import ScriptWriter
    in_file=r"C:\Users\devora\Sources\csm\test_cases\inbal\reading fragments\model-endmdl-withIDS.pdb"
    command_file=r"C:\Users\devora\Sources\temp\csm_multi_test\input\sym.txt"
    out_folder=r"C:\Users\devora\Sources\temp\csm_multi_test\output"

    commands=[]
    dictionary_args_array = []
    with open(command_file, 'r') as file:
        for command in file:
            commands.append(command)
            args = command.split()
            dictionary_args = get_parsed_args(args)
            dictionary_args_array.append(dictionary_args)

    mols=MoleculeReader.multiple_from_file(in_file, use_chains=False)
    results=[]
    for index, mol in enumerate(mols):
        mol_results=[]
        for args in dictionary_args_array:
            args['molecule']=mol
            calc, result, dictionary_args=_run_calculation(args)
            mol_results.append(result)
        results.append((str(index), mol_results))
    format=get_format(None, filename=in_file)
    writer=ScriptWriter(results, format, commands, out_folder)
    writer.write()

def run_no_return(args=[]):
    run(args)


if __name__ == '__main__':
    timer = timeit.Timer(lambda: run(args=sys.argv[1:]))
    time = timer.timeit(number=1)
    print("Runtime:", time, "seconds")
