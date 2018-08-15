import csv

from csm.calculations import Approx, Trivial, Exact, ParallelApprox, DirectionChooser
from csm.calculations.approx.dirs import get_direction_chooser
from csm.calculations.basic_calculations import CalculationTimeoutError
from csm.input_output.readers import read_perm
from csm.input_output.formatters import csm_log as print
from csm.main.normcsm import norm_calc


def do_calculation(command, perms_csv_name=None, parallel=False, print_approx=False, **dictionary_args):
    calc_type=command
    if calc_type=="exact":
        #get perm if it exists:
        dictionary_args['perm'] = read_perm(**dictionary_args)
        csm_state_tracer_func= None
        if perms_csv_name:
            csv_file = open(perms_csv_name, 'w')
            perm_writer = csv.writer(csv_file, lineterminator='\n')
            perm_writer.writerow(['Permutation', 'Direction', 'CSM'])
            csm_state_tracer_func = lambda state: perm_writer.writerow(
                [[p + 1 for p in state.perm],
                 state.dir,
                 state.csm, ])
        calc=Exact(**dictionary_args, callback_func=csm_state_tracer_func)

    if calc_type=="approx":
        dir_chooser = get_direction_chooser(**dictionary_args)
        dictionary_args["direction_chooser"] = dir_chooser
        if parallel:
            calc=ParallelApprox(**dictionary_args)
        else:
            if print_approx:
                def log(self, *args, **kwargs):
                    print(*args)
                dictionary_args["log_func"]=log
            calc = Approx(**dictionary_args)

    if calc_type=="trivial":
        calc = Trivial(**dictionary_args)

    #run the calculation
    calc.calculate(**dictionary_args)
    return calc.result


def single_calculation(molecule, dictionary_args):
    print("Molecule:", molecule.metadata.header())
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
