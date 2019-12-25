import csv
import json
import os
import pytest
import shutil

from csm.main.csm_run import csm_run
from tests.local_settings import test_dir


class RunThings():
    def _run_args(self, args_str, results_folder):
        args_str += " --output {} --overwrite".format(results_folder)
        args = args_str.split()
        tewst=(os.getcwd())
        results_arr = csm_run(args)
        return results_arr


class TestApproxCorrect(RunThings):
    os.chdir(test_dir)
    results_folder = "csm_tests"

    def run_args(self, args_str):
        os.chdir(r'..\..\correctly_tests\files_for_tests')
        return super()._run_args(args_str, self.results_folder)

    def test_1(self):
        file_name = '5nnm.pdb'
        cmd = "approx c2 --input {input} --use-sequence --use-chain".format(input=file_name)
        result = self.run_args(cmd)

        output_extra = open(r"csm_tests\extra.txt").readlines()[1:]
        excepted_extra = open(r"excepted_files\extra_5nnm_test1.txt").readlines()[1:]
        assert output_extra == excepted_extra
        assert result[0][0].csm == pytest.approx(0.0642662295, rel=1e-8)

    def test_2(self):
        file_name = '1aoj.pdb'  # file clean after pdb-prep
        cmd = "approx c2 --input {input} --use-sequence --many-chain".format(input=file_name)
        result = self.run_args(cmd)

        output_extra = open(r"csm_tests\extra.txt").readlines()[1:]
        excepted_extra = open(r"excepted_files\extra_1aoj_test2.txt").readlines()[1:]
        assert output_extra == excepted_extra
        assert result[0][0].chain_perm == [1, 0]
        assert result[0][0].csm == pytest.approx(0.4623219122, rel=1e-8)

    def test_3(self):
        """ the result is ~34, because the flag '--input-chain-perm'. """
        file_name = '1aoj.pdb'  # file clean after pdb-prep
        cmd = "approx c2 --input {input} --use-sequence --input-chain-perm perm_AB_01.txt --use-backbone --print-approx --use-best-dir".format(input=file_name)
        result = self.run_args(cmd)

        output_extra = open(r"csm_tests\extra.txt").readlines()[1:]
        excepted_extra = open(r"excepted_files\extra_1aoj_test3.txt").readlines()[1:]
        count_diff = 0
        for l_exp, l_out in zip(excepted_extra, output_extra):
            if l_exp != l_out: # the lines "Calculating for initial direction:" and the distance.
                count_diff += 1
        assert count_diff <= 6
        assert result[0][0].chain_perm == [0, 1]

    def test_4(self):
        file_name = '1aoj.pdb'  # file clean after pdb-prep
        cmd = "approx c2 --input {input} --select-atoms 17-20,19-21".format(input=file_name)
        result = self.run_args(cmd)

        assert len(result[0][0].molecule._atoms) == 5
        assert result[0][0].csm == pytest.approx(3.450308754, rel=1e-8)

    def test_5(self):
        file_name = '3kbl.pdb'  # chains' lengths: [465, 384, 455, 380]
        cmd = "approx c2 --input {input}  --remove-hy --select-mols 1 --timeout 60 --babel-bond --polar".format(input=file_name)
        result = self.run_args(cmd)

        output_tsv = open(r"csm_tests\approx\3kbl_L01_c2.tsv").readlines()
        excepted_tsv = open(r"excepted_files\3kbl_L01_c2.tsv").readlines()  # 3kbl_L01_c2_test5.tsv

        for excepted_row, output_row in zip(excepted_tsv, output_tsv):
            # remove the column "Runtime"
            excepted_row = excepted_row.split('\t')
            output_row = output_row.split('\t')
            excepted_row = list(excepted_row[i] for i in range(len(excepted_row)) if i != 10)
            output_row = list(output_row[i] for i in range(len(output_row)) if i != 10)
            assert excepted_row == output_row
        assert result[0][0].csm == pytest.approx(1.115601065, rel=1e-8)

    def test_6(self):
        file_name = '2l34.pdb'
        cmd = "approx c2 --input {input} --use-sequence --select-mols 1 --fibonacci 1 --dir 0 0 1 " \
              "--greedy".format(input=file_name)
        result = self.run_args(cmd)

        assert (result[0][0].ongoing_statistics['approx'])['c2'][0]['dir'] == (0,0,1)

    def test_7(self):
        file_name = '6g8o.pdb'
        cmd = "approx c2 --input {input} --use-sequence --select-mols 1 --remove-hy --verbose --connect".format(input=file_name)
        with pytest.raises(ValueError):  # ValueError: --connect work only with .xyz files
            result = self.run_args(cmd)

    def test_8(self):
        file_name = '6g8o.pdb'
        cmd = "approx c2 --input {input} --use-sequence --use-chains --select-mols 1 --remove-hy --verbose  --ignore-atoms 1-10 --ignore-sym".format(input=file_name)
        result = self.run_args(cmd)

        print(len(result[0][0].molecule._atoms))

    def test_9(self):
        file_name = '6g8o.pdb'
        cmd = "approx c2 --input {input} --use-sequence --keep-structure --select-mols 1".format(input=file_name)
        result = self.run_args(cmd)

        assert result[0][0].csm == pytest.approx(0.000033, abs=1e-8)

    def test_10(self):
        file_name = 'for_try/4v2t.pdb'  # your chains' lengths: [135, 474, 135, 135, 474, 135, 135, 474, 135, 135, 475, 135, 135, 474, 135, 135, 474, 135, 135, 474, 135, 135, 474, 135, 135, 474, 135, 135, 474, 135, 135, 474, 135, 135, 474, 135, 474, 135, 135]
        cmd = "approx c13 --input {input} --use-sequence --remove-hy --select-mols 1".format(input=file_name)
        result = self.run_args(cmd)

#--print-approx
"""
    v--connect
    v--select-atoms
    v--ignore-atoms
    v--remove-hy
    v--select-mols
    v--ignore-sym  # ignore symbol
    v--use-mass
    v--babel-bond  # openbabel compute the bonds
    v--use-sequence  # pdb file. work with other file?
    v--use-chains  # approx, trivial
    --read-fragments  # Read fragments from .mol or .pdb file as chains
    v--use-backbone  # other file from pdb?
    
    output:
    --out-format
    v--simple  # Ignores other output flags?
    --legacy-output
    v--overwrite 
    v--verbose  # create one more file. what do file name?
    v--polar # polar coor in statistics created by --verbose
    --json-output # to which file it's print?
    --legacy-files # what mean 'legacy'?
    v--timeout
    --sn-max  # only for chirality  
    
    
    approx:
    ?--detect-outliers  # "Use outlier detection to improve guesses for initial directions in approx algorithm. Only activated for more than 10 equivalence groups."
    --no-orthogonal  # Don't add orthogonal directions to calculated directions
    v--use-best-dir
    v--fibonacci n  # n starting directions
    v--dir  x y z  # specific dir
    
    (approx algorithm:)
    --greedy  # default hungarian
    --many-chains  
    --keep-structure
    --selective  k  # with fibonacci
    --input-chain-perm 
    --print-approx # print the logs to the screen
"""



