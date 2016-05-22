import csv

import os
import glob
import numpy as np

from csm_run import run
from csm.molecule.molecule import Molecule


outputdirectory=r'D:\UserData\devora\Sources\csm\test_cases\inbal\expected_output\expected_output_'
inputdirectory=r'D:\UserData\devora\Sources\csm\test_cases\inbal\input\input_'
testnames= ['1_methane_csm',
            '2_biphenyls',
            '3_cyclopentadiene',
            '4_clusters']



def check_result(result, validate):
    # fetching the result and validation data:
    r_perm = [p + 1 for p in result.perm]
    v_perm = [int(i) for i in validate["perm"].split()]
    r_csm = result.csm
    v_csm = float(validate["csm"])
    r_dir = result.dir
    v_dir = [float(i) for i in validate["dir"].split()]
    r_atoms = result.symmetric_structure
    v_mol = Molecule.from_string(validate["out_atoms"], "xyz", False)
    v_atoms = np.asarray([np.asarray(atom.pos) for atom in v_mol.atoms])

    # build the dictionary:
    res = {}
    res["molecule id"] = validate['mol_index']
    res["symmetry"] = validate['symm']
    res["csm result"] = r_csm
    res["dir result"] = r_dir
    res["perm result"] = r_perm
    res["atoms result"] = r_atoms
    res["csm expected"] = v_csm
    res["dir expected"] = v_dir
    res["perm expected"] = v_perm
    res["atoms expected"] = v_atoms

    # check the status:
    failed = False
    res["status"] = "Failed"
    res["message"] = []

    # check csm:
    if abs(r_csm - v_csm) > .0001:  # greater than 4 decimal points: the validation stops, hence the result will always be bigger (since it continue for several more decimal points
        failed = True
        res["message"].append("CSM mismatch")

    # chck symmetric structure
    atom_pos_not_match = []
    for i in range(len(r_atoms)):
        for index in range(3):
            val1 = r_atoms[i][index]
            val2 = v_atoms[i][index]
            if abs(val1 - val2) > .01:
                atom_pos_not_match.append((i, index, val1, val2))
    if atom_pos_not_match:
        failed = True
        res["message"].append("Structure mismatch")

    res['verify']=[]
    # check perm
    if r_perm != v_perm:
        res['verify']=True
        res["message"].append("Perm mismatch")

    # check dir
    dir_not_match = False
    dir_not_times_minus = False
    for i in range(3):
        val1 = r_dir[i]
        val2 = (v_dir[i])
        if abs(val1 - val2) > .0001:
            dir_not_match = True
            if abs(-val1 - val2) > .0001:
                dir_not_times_minus = True
    if dir_not_match:
        if dir_not_times_minus:
            res["message"].append("Dir mismatch")
        else:
            res["message"].append("Dir negative")

    if failed:
        return atom_pos_not_match, res
    res["status"] = "OK"
    return atom_pos_not_match, res

def parse_and_run(workdir):
        args=[]
        for argfilename in glob.glob(os.path.join(workdir, '*sym*.txt')):
            argfile=open(argfilename, 'r')
            for line in argfile:
                args.append(line)

        observed_results=[]
        for molfilename in glob.glob(os.path.join(workdir, '*.xyz')):
            for arg in args:
                    fixed_args=arg.replace("__IN__", molfilename).replace("__OUT__", "outfile.out").replace("-formatxyz","").replace("-useperm","--useperm")
                    fixed_args=fixed_args.split()
                    result = run(fixed_args[1:])
                    observed_results.append(result)

        return observed_results

def parse_expected(test):
    resf=os.path.join(outputdirectory+test,'csmresults.log')
    resultsfile=open(resf, 'r')
    e = resultsfile.read().split("MOL_")
    e= e[1:]
    results=[]
    for result in e:
        dict={}
        dict['mol_index']=''


def run_tests():
    for test in testnames:
        workdir=inputdirectory+test
        os.chdir(workdir)

        observed_results=''#parse_and_run(workdir)
        expected_results=parse_expected(test)


        with open(test+'.csv', 'w') as csvfile:
            fieldnames = ['molecule id', 'symmetry', 'status', 'message', 'verify', 'csm result', 'csm expected', 'atoms result',
                      'atoms expected', 'perm result', 'perm expected', 'dir result', 'dir expected']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, lineterminator='\n')
            writer.writeheader()
            print(len(observed_results))
            print("----")
            for i in range(len(observed_results)):
                observed=observed_results[i]
                expected=expected_results[i]
                row=check_result(observed, expected)
                writer.writerow(row)








run_tests()