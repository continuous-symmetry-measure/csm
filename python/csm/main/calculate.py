import csv

from csm.calculations import Approx, Trivial, Exact, ParallelApprox, DirectionChooser
from csm.calculations.approx.dirs import get_direction_chooser
from csm.input_output.readers import read_perm


def do_calculation(**dictionary_args):
    calc_type=dictionary_args["command"]
    if calc_type=="exact":
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

    if calc_type=="approx":
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

    if calc_type=="trivial":
        calc = Trivial(**dictionary_args)

    #run the calculation
    calc.calculate(**dictionary_args)
    return calc.result