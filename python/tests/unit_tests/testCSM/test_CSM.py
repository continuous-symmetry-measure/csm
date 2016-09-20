import pytest
import json
from csm.main.csm_run import run
from tests.utils.run_test import close_enough


def run_test(args, molecule):
    args[1]=molecule
    args[2]=r'D:\UserData\devora\Sources\csm\python\tests\unit_tests\noonecares.txt'
    result=run(args)
    return result


def test_csm(args, molecule, expected, equiv):
    result=run_test(args, molecule)
    assert close_enough(result.csm, expected.csm)

def test_symmetric_structure(args, molecule, expected, equiv):
    result=run_test(args, molecule)
    atom_pos_not_match = []
    for i in range(len(expected.molecule)):
        for index in range(3):
            val1 = expected.symmetric_structure[i][index]
            val2 = result.symmetric_structure[i][index]
            if not close_enough(val1, val2):
                atom_pos_not_match.append((str(i), index, val1, val2))
    assert(len(atom_pos_not_match)==0)

def test_perm(args, molecule, expected, equiv):
    result=run_test(args, molecule)
    symm=args[0]

    exp_p=expected.perm
    res_p=result.perm
    if expected.perm != result.perm:
        if symm in equiv:
            for pair in equiv[symm]:
                if expected.perm in pair and result.perm in pair:
                    res_p=exp_p
                    break

    assert(exp_p==res_p)
