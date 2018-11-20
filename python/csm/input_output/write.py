import json
import sys

from csm.calculations.data_classes import CSMResult
from csm.input_output.formatters import csm_log as print
from csm.input_output.formatters import format_CSM
from csm.input_output.readers import read_from_sys_std_in
from csm.input_output.writers import ScriptWriter, LegacyFormatWriter


def write_results(results, in_format=None, out_format=None, simple=False, legacy=False, pipe=False, **kwargs):
    '''
    :param results: must be an array of arrays, outer array molecules, inner array commands
    :param kwargs:
    :return:
    '''
    results_arr = results
    if simple:
        for mol_index, mol_result in enumerate(results_arr):
            for lin_index, line_result in enumerate(mol_result):
                print("mol", line_result.molecule.metadata.appellation(), "cmd", lin_index + 1, " CSM: ",
                      format_CSM(line_result.csm))
        return

    if pipe:
        sys.stdout.write(
            json.dumps([[result.to_dict() for result in mol_results_arr] for mol_results_arr in results_arr], indent=4))
        return

    if out_format:
        format = out_format
    elif in_format:
        format = in_format
    else:
        format = results_arr[0][0].molecule.metadata.format

    if legacy:
        if len(results_arr) > 1 or len(results_arr[0]) > 1:
            raise ValueError("Legacy result writing only works for a single molecule and single command")
        result = results_arr[0][0]
        writer = LegacyFormatWriter(result, format)
        with open(kwargs["out_file_name"], 'w') as f:
            writer.write(f)
        return

    writer = ScriptWriter(results_arr, format, **kwargs)
    writer.write()


def write(**dictionary_args):
    raw_json = read_from_sys_std_in()
    less_raw_json = json.loads(raw_json)
    results = [[CSMResult.from_dict(result_dict) for result_dict in mol_arr] for mol_arr in less_raw_json]
    write_results(results, **dictionary_args)
    return
