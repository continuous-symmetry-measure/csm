import os
import pytest
from csm.main.csm_run import csm_run as csmrun
from tests.utils.run_test import close_enough, YaffaError
from conftest import test_folder, output_file, my_tolerance


os.chdir(test_folder)

@pytest.mark.parametrize("args, expected",[
    (['c2', r'--input test1\AgCu10p1.xyz', "--output " +output_file], 6.5620),
    (['cs', r'--input test2\ZnCu10-p1.csm', "--output " +output_file], 3.3066),
    (['c4', r'--input test3\c_in_1282148276_benzene.mol', "--output " +output_file], 41.6651),
    (['c2', r'--input test4\ZnCu10-p9.csm', "--output " +output_file, '--ignore-sym'], 0.0200),
    (['s4', r'--input test6\AuCu10p6.xyz', "--output " +output_file, '--babel-bond'], 100.0000),
    #(['s8', 'test7\AuCu10p9.xyz', "--output " +output_file], 24.0465),
    (['c3', r'--input test8\h2o2.xyz', "--output " +output_file], 7.1183),
    (['c3', r'--input test9\h2o2.xyz', "--output " +output_file, '--remove-hy'], 0.00),
    #(['ci', r'--input test10\PHI.cif', "--output " +output_file], 6.5620),
    #(['cs', r'--input test11\PHI.cif', "--output " +output_file], 38.8953),
    #(['ci', r'--input test12\CUBLEQ.CIF', "--output " +output_file, '--remove-hy', '--babel-bond'], "NO RESULT IN FOLDER"),
    #(['ci', r'--input test13\CUBLEQ.PDB', "--output " +output_file, '--remove-hy', '--babel-bond'], 33.2130),
    #(['ci', r'--input test13csm\CUBLEQ.csm', "--output " +output_file, '--remove-hy', '--babel-bond'], 33.2130),
    #(['ci', r'--input test14\AFADOA.PDB', "--output " +output_file, '--remove-hy', '--babel-bond'], 4.6905),
    #(['ch', r'--input test15\AGA-alpha.pdb', "--output " +output_file],  8.8716),
    (['c4', r'--input test16\tetralin.mol', "--output " +output_file], 20.7659),
    (['s6', r'--input test17\alpha-D-glucopyranose.mol', "--output " +output_file], 98.9854),
    #(['c4', r'--input test18\phthalocyanin.mol', "--output " +output_file, '--remove-hy'], 29.0559),
    #(['c2', r'--input test19\11Bg.pdb', "--output " +output_file, '--approx', '--babel-bond'], 0.1505),
    (['ch', r'--input test20\h2o2-twist1.xyz', "--output " +output_file, '--babel-bond'], 2.3658),
    (['ch', r'--input test21\h2o2-twist2.xyz', "--output " +output_file, '--babel-bond'], 4.3314),
    (['ch', r'--input test22\h2o2-twist3.xyz', "--output " +output_file, '--babel-bond'], 6.2427),
    (['ch', r'--input test23\bis(dth)copper(I).mol', "--output " +output_file], 0.000),
    (['ch', r'--input test24\bis(dth)copper(I)-dist.mol', "--output " +output_file], 3.0302),

                                            ])
def test_exact(args, expected):
        args=["exact"]+args
        print(args)
        result = csmrun(args)
        if not close_enough(result.csm, result.formula_csm):
            print("Yaffa!")
            raise YaffaError
        assert close_enough(result.csm, expected)



