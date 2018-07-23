import json
import sys

import os

from csm.calculations.data_classes import CSMResult
from csm.input_output.formatters import format_CSM
from csm.input_output.readers import read_from_sys_std_in
from csm.input_output.writers import OldFormatFileWriter, ScriptWriter
from csm.input_output.formatters import csm_log as print

def write_results(results, **kwargs):
    '''
    :param results: must be an array of arrays, outer array molecules, inner array commands
    :param kwargs:
    :return:
    '''
    results_arr=results
    if kwargs['simple']:
        for mol_index, mol_result in enumerate(results_arr):
            for lin_index, line_result in enumerate(mol_result):
                print("mol", mol_index+1, "cmd", lin_index+1, " CSM: ", format_CSM(line_result.csm))
        return

    if kwargs["pipe"]:
        sys.stdout.write(json.dumps([[result.to_dict() for result in mol_results_arr] for mol_results_arr in results_arr], indent=4))
        return

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
        format=results_arr[0][0].molecule.metadata.format
    writer = ScriptWriter(results_arr, format, **kwargs)
    writer.write()





def write(**dictionary_args):
    raw_json = read_from_sys_std_in()
    less_raw_json = json.loads(raw_json)
    results = [[CSMResult.from_dict(result_dict) for result_dict in mol_arr] for mol_arr in less_raw_json]
    write_results(results, **dictionary_args)
    return