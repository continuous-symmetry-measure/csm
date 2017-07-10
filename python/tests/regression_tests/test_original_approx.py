

import os
import pytest
from csm.main.csm_run import run as csmrun
from tests.utils.run_test import close_enough, YaffaError

output_file=r'C:\Users\devora.CHELEM\Sources\temp\csm_tests_output.txt'
testdir=r'C:\Users\devora.CHELEM\Sources\csm\test_cases\old_test_cases\original_test_cases'
os.chdir(testdir)

@pytest.mark.parametrize("args, expected",[
    (['c2', r'test1\AgCu10p1.xyz', output_file], 6.5620),
    (['cs', r'test2\ZnCu10-p1.csm', output_file], 3.3066),
    (['c4', r'test3\c_in_1282148276_benzene.mol', output_file], 41.6651),
    (['c2', r'test4\ZnCu10-p9.csm', output_file, '--ignore-sym'], 0.0200),
    (['s4', r'test6\AuCu10p6.xyz', output_file, '--babel-bond'], 100.0000),
    (['s8', 'test7\AuCu10p9.xyz', output_file], 24.0465),
    (['c3', r'test8\h2o2.xyz', output_file], 7.1183),
    (['c3', r'test9\h2o2.xyz', output_file, '--remove-hy'], 0.00),
    (['ci', r'test10\PHI.cif', output_file], 6.5620),
    (['cs', r'test11\PHI.cif', output_file], 38.8953),
    #(['ci', r'test12\CUBLEQ.CIF', output_file, '--remove-hy', '--babel-bond'], "NO RESULT IN FOLDER"),
    (['ci', r'test13\CUBLEQ.PDB', output_file, '--remove-hy', '--babel-bond'], 33.2130),
    (['ci', r'test13csm\CUBLEQ.csm', output_file, '--remove-hy', '--babel-bond'], 33.2130),
    (['ci', r'test14\AFADOA.PDB', output_file, '--remove-hy', '--babel-bond'], 4.6905),
    (['ch', r'test15\AGA-alpha.pdb', output_file],  8.8716),
    (['c4', r'test16\tetralin.mol', output_file], 20.7659),
    (['s6', r'test17\alpha-D-glucopyranose.mol', output_file], 98.9854),
    (['c4', r'test18\phthalocyanin.mol', output_file, '--remove-hy'], 29.0559),
    (['c2', r'test19\11Bg.pdb', output_file, '--babel-bond'], 0.1505),
    (['ch', r'test20\h2o2-twist1.xyz', output_file, '--babel-bond'], 2.3658),
    (['ch', r'test21\h2o2-twist2.xyz', output_file, '--babel-bond'], 4.3314),
    (['ch', r'test22\h2o2-twist3.xyz', output_file, '--babel-bond'], 6.2427),
    (['ch', r'test23\bis(dth)copper(I).mol', output_file], 0.000),
    (['ch', r'test24\bis(dth)copper(I)-dist.mol', output_file], 3.0302),
                                            ])
def xtest_approx(args, expected):
        print(args)
        print(expected)
        args.append('--approx')
        #args.append('--no-hungarian')
        #args.append('--new-chains')
        result = csmrun(args)
        if not close_enough(result.csm, result.formula_csm):
            print("Yaffa!")
            raise YaffaError
        assert close_enough(result.csm, expected)