'''
A set of tests to check that a variety of input arguments are still valid and haven't broken from changes introduced to the code-- doesn't test validity of result, just
that the code runs successfully and returns A result.
'''


import pytest
from csm.main.csm_run import run

class Runner:
    def run_args(self, args_str):
        args=args_str.split()
        results_arr=run(args)
        result=results_arr[0][0]
        assert result.csm is not None

class TestExact(Runner):
    def test_plain(self):
        self.run_args(r"exact c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt")

    @pytest.mark.parametrize("run_str",[
        r"exact ch --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --sn-max 6",
        r"exact c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2RLA-3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --remove-hy",
        r"exact c3 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --use-mass",
        r"exact s8 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --ignore-sym",
        r"exact s6 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --read-fragments",
        r"exact s2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2RLA-3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --use-sequence",
        r"exact ci --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --use-chains",
        r"exact c4 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --babel-bond",
    ])
    def test_input_args(self, run_str):
        self.run_args(run_str)
    def test_exact_args(self):
        self.run_args(r"exact c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --keep-structure")

    @pytest.mark.parametrize("run_str",
                             [
        r"exact c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --print-denorm --legacy",
        r"exact c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --output-branches --legacy",
        r"exact c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --output-perms perms.csv --legacy",
        r"exact c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --print-local --legacy",
        r"exact c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --json-output --legacy",
        r"exact c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --simple",
    ])
    def test_output_args(self, run_str):
        self.run_args(run_str)

class TestApprox(Runner):
    def test_plain(self):
        self.run_args(r"approx c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt ")
        self.run_args(r"approx c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt  --print-approx")
    @pytest.mark.parametrize("run_str",
                             [

                                 r"approx ch --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt  --sn-max 6",
                                 r"approx ci --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt  --use-chains",
                                 r"approx s2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt  --use-sequence",
                                 r"approx s6 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt  --read-fragments",

                             ])
    def test_input_args(self, run_str):
        self.run_args(run_str) 
    @pytest.mark.parametrize("run_str",
                             [
                                 r"approx c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\approx\csm_testoutput.txt  --detect-outliers",
                                 r"approx c3 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\approx\csm_testoutput.txt  --no-orthogonal",
                                 r"approx c5 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\approx\csm_testoutput.txt  --detect-outliers --no-orthogonal",
                                 r"approx cs --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\approx\csm_testoutput.txt  --fibonacci 20",
                                 r"approx c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\approx\csm_testoutput.txt  --use-best-dir",
                                 r"approx c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\approx\csm_testoutput.txt  --many-chains",
                                 r"approx c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\approx\csm_testoutput.txt  --greedy",
                                 r"approx c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\approx\csm_testoutput.txt  --fibonacci 20 --selective 3",
                                 r"approx c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\approx\csm_testoutput.txt  --statistics C:\Users\devora\Sources\temp\csm_argument_tests\approx\csm_testoutputstats.txt",
                                 r"approx c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\approx\csm_testoutput.txt  --statistics C:\Users\devora\Sources\temp\csm_argument_tests\approx\csm_testoutputstats.txt --polar",
                             ])
    def test_approx_args(self, run_str):
        self.run_args(run_str)
    
    @pytest.mark.parametrize("run_str",
                             [
                                 r"approx c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\approx --output C:\Users\devora\Sources\temp\csm_argument_tests\approx\multi  --use-best-dir",
                                 r"approx c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\approx --output C:\Users\devora\Sources\temp\csm_argument_tests\approx\multi  --fibonacci 3",
                             ])
    def test_multi_mols(self, run_str):
        self.run_args(run_str)
        
    def test_output_args(self):
        self.run_args(r"approx c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt  --json-output --legacy")

class TestTrivial(Runner):
    def test_plain(self):
        self.run_args(r"trivial c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt ")
    
    @pytest.mark.parametrize("run_str", [

        r"trivial ch --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt  --sn-max 6",
        r"trivial c4 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt  --babel-bond",
        r"trivial ci --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt  --use-chains",
        r"trivial s2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt  --use-sequence",
        r"trivial s6 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt  --read-fragments",
    ])
    def test_input_args(self, run_str):
        self.run_args(run_str)

    @pytest.mark.parametrize("run_str", [

        r"trivial ch --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests --output C:\Users\devora\Sources\temp\csm_argument_tests\multi_molecules  --sn-max 6",
        r"trivial c4 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests --output C:\Users\devora\Sources\temp\csm_argument_tests\multi_molecules  --babel-bond",
        r"trivial ci --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests --output C:\Users\devora\Sources\temp\csm_argument_tests\multi_molecules  --use-chains",
        r"trivial s2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests --output C:\Users\devora\Sources\temp\csm_argument_tests\multi_molecules  --use-sequence",
        r"trivial s6 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests --output C:\Users\devora\Sources\temp\csm_argument_tests\multi_molecules  --read-fragments",
    ])
    def test_multiple_mols(self, run_str):
        self.run_args(run_str)