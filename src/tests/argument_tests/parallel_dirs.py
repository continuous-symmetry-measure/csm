import os
from csm.main.csm_run import csm_run


def parallel_dirs_test_indepentent():
    args = ["approx", "c2", "--input",  os.path.join(os.getcwd(), "tests", "argument_tests", "files_for_tests", "3alb-gkt4-h.pdb"), "--use-sequence", "--parallel-dirs"]
    results_arr = csm_run(args)
    if round(results_arr[0][0].csm, 6) != 0.131189:
        raise Exception("Parallel dir test failed!")

if __name__ == "__main__":
    parallel_dirs_test_indepentent()