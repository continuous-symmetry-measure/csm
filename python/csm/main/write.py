import json
import sys

import os

from csm.calculations.data_classes import CSMResult
from csm.input_output.readers import read_from_sys_std_in
from csm.input_output.writers import OldFormatFileWriter, ScriptWriter


def write_results(results, **kwargs):
    try:
        if not isinstance(results[0], list):  # results is a single array
            results = [results]
    except TypeError:  # results isn't an array at all
        results = [[results]]
    except IndexError: #results is an array, and the array is empty
        raise ValueError("Can't write empty results")

    results_arr=results
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
            #if len(results_arr) == 1 and len(results_arr[0]) == 1:
            #    print("You are running a single file and command. Did you want to print to the old format, --legacy?")
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


def write(**dictionary_args):
    raw_json = read_from_sys_std_in()
    less_raw_json = json.loads(raw_json)
    results = [[CSMResult.from_dict(result_dict) for result_dict in mol_arr] for mol_arr in less_raw_json]
    write_results(results, **dictionary_args)
    return