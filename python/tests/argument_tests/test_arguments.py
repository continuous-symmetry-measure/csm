import pytest
from csm.main.csm_run import run

class Runner:
    def run_args(self, args_str):
        args=args_str.split()
        result=run(args)



class TestExact(Runner):
    def test_plain(self):
        self.run_args(r"c2 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol C:\Users\devora\Sources\temp\csm_testoutput.txt")

    @pytest.mark.parametrize("run_str",[
        r"ch C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol C:\Users\devora\Sources\temp\csm_testoutput.txt --sn-max 6",
        r"c2 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol C:\Users\devora\Sources\temp\csm_testoutput.txt --remove-hy",
        r"c3 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol C:\Users\devora\Sources\temp\csm_testoutput.txt --use-mass",
        r"s8 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol C:\Users\devora\Sources\temp\csm_testoutput.txt --ignore-sym",
        r"s6 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol C:\Users\devora\Sources\temp\csm_testoutput.txt --read-fragments",
        r"s2 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2RLA-3.pdb C:\Users\devora\Sources\temp\csm_testoutput.txt --use-sequence",
        r"ci C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol C:\Users\devora\Sources\temp\csm_testoutput.txt --use-chains",
        r"c4 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol C:\Users\devora\Sources\temp\csm_testoutput.txt --babel-bond",
    ])
    def test_input_args(self, run_str):
        self.run_args(run_str)
    def test_exact_args(self):
        self.run_args(r"c2 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol C:\Users\devora\Sources\temp\csm_testoutput.txt --keep-structure")
        self.run_args(r"c2 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol C:\Users\devora\Sources\temp\csm_testoutput.txt --no-constraint")

    @pytest.mark.parametrize("run_str",
                             [
        r"c2 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol C:\Users\devora\Sources\temp\csm_testoutput.txt --print-denorm",
        r"c2 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol C:\Users\devora\Sources\temp\csm_testoutput.txt --output-branches",
        r"c2 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol C:\Users\devora\Sources\temp\csm_testoutput.txt --output-perms perms.csv",
        r"c2 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol C:\Users\devora\Sources\temp\csm_testoutput.txt --print-local",
        r"c2 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol C:\Users\devora\Sources\temp\csm_testoutput.txt --json-output",
    ])
    def test_output_args(self, run_str):
        self.run_args(run_str)

class TestApprox(Runner):
    def test_plain(self):
        self.run_args(r"c2 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb C:\Users\devora\Sources\temp\csm_testoutput.txt --approx")
        self.run_args(r"c2 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb C:\Users\devora\Sources\temp\csm_testoutput.txt --approx --print-approx")
    @pytest.mark.parametrize("run_str",
                             [

                                 r"ch C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb C:\Users\devora\Sources\temp\csm_testoutput.txt --approx --sn-max 6",
                                 r"ci C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb C:\Users\devora\Sources\temp\csm_testoutput.txt --approx --use-chains",
                                 r"s2 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb C:\Users\devora\Sources\temp\csm_testoutput.txt --approx --use-sequence",
                                 r"s6 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb C:\Users\devora\Sources\temp\csm_testoutput.txt --approx --read-fragments",

                             ])
    def test_input_args(self, run_str):
        self.run_args(run_str) 
    @pytest.mark.parametrize("run_str",
                             [
                                 r"c2 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb C:\Users\devora\Sources\temp\csm_testoutput.txt --approx --detect-outliers",
                                 r"c3 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb C:\Users\devora\Sources\temp\csm_testoutput.txt --approx --no-orthogonal",
                                 r"c5 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb C:\Users\devora\Sources\temp\csm_testoutput.txt --approx --detect-outliers --no-orthogonal",
                                 r"cs C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb C:\Users\devora\Sources\temp\csm_testoutput.txt --approx --fibonacci 20",
                                 r"c2 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb C:\Users\devora\Sources\temp\csm_testoutput.txt --approx --use-best-dir",
                                 r"c2 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb C:\Users\devora\Sources\temp\csm_testoutput.txt --approx --many-chains",
                                 r"c2 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb C:\Users\devora\Sources\temp\csm_testoutput.txt --approx --greedy",
                                 r"c2 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb C:\Users\devora\Sources\temp\csm_testoutput.txt --approx --fibonacci 20 --selective 3",
                                 r"c2 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb C:\Users\devora\Sources\temp\csm_testoutput.txt --approx --statistics C:\Users\devora\Sources\temp\csm_testoutputstats.txt",
                                 r"c2 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb C:\Users\devora\Sources\temp\csm_testoutput.txt --approx --statistics C:\Users\devora\Sources\temp\csm_testoutputstats.txt --polar",

                             ])
    def test_approx_args(self, run_str):
        self.run_args(run_str)
        
    def test_output_args(self):
        self.run_args(r"c2 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb C:\Users\devora\Sources\temp\csm_testoutput.txt --approx --json-output")

class TestTrivial(Runner):
    def test_plain(self):
        self.run_args(r"c2 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol C:\Users\devora\Sources\temp\csm_testoutput.txt --trivial")

    def test_input_args(self):
        self.run_args(r"ch C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol C:\Users\devora\Sources\temp\csm_testoutput.txt --trivial --sn-max 6")
        self.run_args(r"c4 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol C:\Users\devora\Sources\temp\csm_testoutput.txt --trivial --babel-bond")
        self.run_args(r"ci C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol C:\Users\devora\Sources\temp\csm_testoutput.txt --trivial --use-chains")
        self.run_args(r"s2 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol C:\Users\devora\Sources\temp\csm_testoutput.txt --trivial --use-sequence")
        self.run_args(r"s6 C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol C:\Users\devora\Sources\temp\csm_testoutput.txt --trivial --read-fragments")