'''
test fo any bugs found, to make sure they don't recur
'''
import os
import shutil
from csm.main.csm_run import csm_run



class RunThings():
    def _run_args(self, args_str, results_folder):
        os.chdir(r'..\..\regression_tests\files_for_tests')
        args_str += " --output {} --overwrite".format(results_folder)
        args = args_str.split()
        results_arr = csm_run(args)
        return results_arr


class Test_0_22_3(RunThings):
    test_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "files_for_tests")
    os.chdir(test_dir)
    results_folder = "csm_tests"
    try:
        shutil.rmtree(results_folder)
    except FileNotFoundError:
        pass
    os.mkdir(results_folder)

    def run_args(self, args_str):
        return super()._run_args(args_str, self.results_folder)

    def test_csm(self):
        #csm files crashed because result writing assumed existence of OBM
        cmd="exact c5 --input nobonds.csm"
        self.run_args(cmd)
    