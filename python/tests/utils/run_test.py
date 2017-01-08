import sys
import json
import math
import os
import warnings

from csm.main.csm_run import run
from csm.molecule.molecule import Molecule

class TestFailedException(Exception):
    pass

class YaffaError(Exception):
    pass

class Expected:
    def __init__(self, dict):
        self.csm=dict['csm']
        self.dir=dict['dir']
        self.perm=[x-1 for x in dict['perm']]
        self.symmetric_structure=Molecule.from_string(dict['symmetric'], 'xyz', initialize=False)._Q
        self.molecule=Molecule.from_string(dict['normalized'], 'xyz')

def close_enough(val1, val2):
    test=math.fabs(val1-val2)
    if test>.0001:
        return False
    return True




def check(result, exp, equiv, test_folder_output, symm):
    try:
        expected=Expected(exp)

        #csm:
        if not close_enough(expected.csm, result.csm):
            raise TestFailedException("Mismatching CSMs")

        #perm:
        mystr = "IDENTICAL"
        if expected.perm != result.perm:
            mystr="MISMATCH"
            if symm in equiv:
                for pair in equiv[symm]:
                    if expected.perm in pair and result.perm in pair:
                        mystr = "EQUIVALENT"

        #symmetric structure:
        atom_pos_not_match = []
        for i in range(len(expected.molecule)):
            for index in range(3):
                val1 = expected.symmetric_structure[i][index]
                val2 = result.symmetric_structure[i][index]
                if not close_enough(val1, val2):
                    atom_pos_not_match.append((str(i), index, val1, val2))
        if atom_pos_not_match:
            warnings.warn("Mismatching symmetric structure")

        if not close_enough(result.formula_csm, result.csm):
            raise YaffaError("Yaffa CSM: %f, calculated CSM: %d" % (result.formula_csm, result.csm))

    finally:
            with open(os.path.join(test_folder_output, 'passed'), 'a') as f:
                f.write("\n\t\tcsm:")
                f.write("\n\t\tperm: "+mystr+str([x+1 for x in expected.perm])+ " "+str([x+1 for x in result.perm]))
                f.write("\tdir: "+str(expected.dir)+ " " + str(result.dir))
                f.write("\n\t\tsymmetric structure mismatches: " + str(atom_pos_not_match))




def run_test(test_folder):
    with open(os.path.join(test_folder, 'input.json')) as f:
        input_dict=json.load(f)
    with open(os.path.join(test_folder, 'output.json')) as f:
        output_dict = json.load(f)
    molecule_folder = os.path.join(test_folder, 'molecules')
    if not os.path.isdir(molecule_folder):
        warnings.warn("No molecule folder for test %s" % test_folder)
        return

    test_output_folder = os.path.join(test_folder, 'results')
    os.makedirs(test_output_folder, exist_ok=True)
    with open(os.path.join(test_output_folder, 'passed.out'), 'w') as f:
        f.write(test_folder)

    for key in input_dict['runs']:
        args=input_dict['runs'][key]
        args[2] = os.devnull
        with open(os.path.join(test_output_folder, 'passed.out'), 'a') as f:
            f.write("\ntest: " + str(key))

        for molecule in os.listdir(molecule_folder):
            try:
                args[1]=os.path.join(molecule_folder, molecule)
                result=run(args)
                mol_index=molecule.split(".")[0]
                with open(os.path.join(test_output_folder, 'passed.out'), 'a') as f:
                    f.write("\n\tmolecule:" + str(mol_index))
                expected=output_dict[key][mol_index]
                check(result, expected, input_dict['equiv_perms'], test_folder_output, args[0])

            except TestFailedException:
                print("FAILED:", key, mol_index)
                raise

            except YaffaError:
                with open(os.path.join(test_output_folder, 'yaffa.out'), 'a') as f:
                    f.write("\n\n\nYAFFA:"+ str(key)+ " "+ str(mol_index)+ " " + str(result.csm)+ " " + str(result.yaffa_csm))










if __name__ == '__main__':
    args = sys.argv[1:]
    test_folder=args[0]
    run_test(test_folder)