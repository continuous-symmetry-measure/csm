import sys
import json

import math

import os
from csm.main.csm_run import run
from csm.molecule.molecule import Molecule

class TestFailedException(Exception):
    pass

class Expected:
    def __init__(self, dict):
        self.csm=dict['csm']
        self.dir=dict['dir']
        self.perm=dict['perm']
        self.symmetric_structure=Molecule.from_string(dict['symmetric'], 'xyz')._Q
        self.molecule=Molecule.from_string(dict['normalized'], 'xyz')

def close_enough(val1, val2):
    test=math.fabs(val1-val2)
    if test>.0001:
        return False
    return True

def check(result, exp):
    expected=Expected(exp)
    if close_enough(expected.csm, result.csm):
        pass
    else:
        raise TestFailedException


def run_test(test_folder):
    print(os.getcwd())
    with open(os.path.join(test_folder, 'input.json')) as f:
        input_dict=json.load(f)
    with open(os.path.join(test_folder, 'output.json')) as f:
        output_dict = json.load(f)

    molecule_folder =os.path.join(test_folder, 'molecules')

    for key in input_dict['runs']:
        args=input_dict['runs'][key]
        args[2]=os.path.join(r'D:\UserData\devora\Sources\csm\python\tests','whocares.out')

        for molecule in os.listdir(molecule_folder):
            try:
                args[1]=os.path.join(molecule_folder, molecule)
                result=run(args)
                mol_index=molecule.split(".")[0]
                expected=output_dict[key][mol_index]
                check(result, expected)
                print("PASSED")
            except TestFailedException:
                print("FAILED:", key, mol_index)
                raise

    with open(os.path.join(test_folder, 'passed'), 'w') as f:
        f.write("passed")







if __name__ == '__main__':
    args = sys.argv[1:]
    test_folder=args[0]
    run_test(test_folder)