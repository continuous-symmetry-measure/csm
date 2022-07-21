'''
A set of tests to check that a variety of input arguments are still valid and haven't broken from changes introduced to the code-- doesn't test validity of result, just
that the code runs successfully and returns A result.
'''
import csv
import json

import os
import pytest
import shutil

from csm.main.csm_run import csm_run, calc, get_parsed_args
from tests.test_settings import test_dir

test_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "files_for_tests")

class RunThings():
    def _run_args(self, args_str, results_folder):
        args_str += " --output {} --overwrite".format(results_folder)
        args = args_str.split()
        tewst=(os.getcwd())
        results_arr = csm_run(args)
        return results_arr


class TestBasic(RunThings):
    os.chdir(test_dir)
    results_folder = "csm_tests"

    # try:
    #    shutil.rmtree(results_folder)
    # except FileNotFoundError:
    #    pass
    # os.mkdir(results_folder)

    def run_args(self, args_str):
        curr_file_dir = os.path.dirname(os.path.realpath(__file__))
        os.chdir(os.path.join(curr_file_dir, 'files_for_tests'))
        return super()._run_args(args_str, self.results_folder)

    # input
    def test_connect(self):
        # --connect reads xyz connectivity file, default connectivity.txt

        # baseline:
        # cmd = "exact c2 --input ferrocene.xyz"
        # results=self.run_args(cmd)
        # assert len(results[0][0].molecule._bondset)==0

        # default:
        cmd = "exact c2 --input ferrocene.xyz --connect"
        results = self.run_args(cmd)
        assert len(results[0][0].molecule.bondset) == 40

        # specified:
        cmd = "exact c2 --input ferrocene.xyz --connect connectivitycopy.txt"
        results = self.run_args(cmd)
        assert len(results[0][0].molecule.bondset) == 40

    def test_remove_hy(self):

        # baseline:
        # cmd="exact cs --input 4-helicene.mol --keep-structure"
        # results=self.run_args(cmd)
        # assert len(results[0][0].molecule) == 30

        cmd = "exact cs --input 4-helicene.mol --remove-hy"
        results = self.run_args(cmd)
        assert len(results[0][0].molecule) == 18

        # test output
        with open(os.path.join(self.results_folder, "resulting_symmetric_coordinates.mol"), 'r') as file:
            file.readline()
            file.readline()
            file.readline()# skip the first 3 lines, which includes the openbabel line, which changes every time
            output_str = file.read()

        expected_openbabel24='''18 21  0  0  0  0  0  0  0  0999 V2000
   -0.9393   -0.9701    0.0521 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2312   -0.4645    0.0947 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3877   -1.2518    0.1527 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2977   -2.6016    0.1706 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0472   -3.1637    0.1304 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8825   -2.3782    0.0720 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3456   -3.0166    0.0333 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5227   -2.2960   -0.0247 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5077   -0.9072   -0.0456 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2892   -0.1740   -0.0088 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4152    1.2854   -0.0363 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6842    1.8958   -0.0955 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8236    1.1097   -0.1285 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7400   -0.2690   -0.1041 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8308    3.2939   -0.1227 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7482    4.1356   -0.0935 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4848    3.5819   -0.0368 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6370    2.1903   -0.0093 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  6  2  0  0  0  0
  1 10  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  6  7  1  0  0  0  0
  7  8  2  0  0  0  0
  8  9  1  0  0  0  0
  9 10  1  0  0  0  0
  9 14  2  0  0  0  0
 10 11  2  0  0  0  0
 11 12  1  0  0  0  0
 11 18  1  0  0  0  0
 12 13  2  0  0  0  0
 12 15  1  0  0  0  0
 13 14  1  0  0  0  0
 15 16  2  0  0  0  0
 16 17  1  0  0  0  0
 17 18  2  0  0  0  0
M  END

$$$$
'''
        expected_openbabel30='''
 18 21  0  0  0  0  0  0  0  0999 V2000
   -0.9393   -0.9701    0.0521 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2312   -0.4645    0.0947 C   0  0  0  0  0  3  0  0  0  0  0  0
   -3.3877   -1.2518    0.1527 C   0  0  0  0  0  3  0  0  0  0  0  0
   -3.2977   -2.6016    0.1706 C   0  0  0  0  0  3  0  0  0  0  0  0
   -2.0472   -3.1637    0.1304 C   0  0  0  0  0  3  0  0  0  0  0  0
   -0.8825   -2.3782    0.0720 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3456   -3.0166    0.0333 C   0  0  0  0  0  3  0  0  0  0  0  0
    1.5227   -2.2960   -0.0247 C   0  0  0  0  0  3  0  0  0  0  0  0
    1.5077   -0.9072   -0.0456 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2892   -0.1740   -0.0088 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4152    1.2854   -0.0363 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6842    1.8958   -0.0955 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8236    1.1097   -0.1285 C   0  0  0  0  0  3  0  0  0  0  0  0
    2.7400   -0.2690   -0.1041 C   0  0  0  0  0  3  0  0  0  0  0  0
    1.8308    3.2939   -0.1227 C   0  0  0  0  0  3  0  0  0  0  0  0
    0.7482    4.1356   -0.0935 C   0  0  0  0  0  3  0  0  0  0  0  0
   -0.4848    3.5819   -0.0368 C   0  0  0  0  0  3  0  0  0  0  0  0
   -0.6370    2.1903   -0.0093 C   0  0  0  0  0  3  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  6  2  0  0  0  0
  1 10  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  6  7  1  0  0  0  0
  7  8  2  0  0  0  0
  8  9  1  0  0  0  0
  9 10  1  0  0  0  0
  9 14  2  0  0  0  0
 10 11  2  0  0  0  0
 11 12  1  0  0  0  0
 11 18  1  0  0  0  0
 12 13  2  0  0  0  0
 12 15  1  0  0  0  0
 13 14  1  0  0  0  0
 15 16  2  0  0  0  0
 16 17  1  0  0  0  0
 17 18  2  0  0  0  0
M  END

$$$$

        '''
        assert output_str.strip()==expected_openbabel24.strip() or output_str.strip()==expected_openbabel30.strip()

    def test_select_atoms(self):
        # --select-atoms removes specific atoms.
        # baseline:
        # cmd="exact c2 --input 4-helicene.mol --keep-structure"
        # results=self.run_args(cmd)
        # assert len(results[0][0].molecule) == 30

        cmd = "exact c2 --input 4-helicene.mol --select-atoms 15-19,1,2"
        results = self.run_args(cmd)
        assert len(results[0][0].molecule) == 7

        # test output
        with open(os.path.join(self.results_folder, "resulting_symmetric_coordinates.mol"), 'r') as file:
            file.readline()
            file.readline()
            file.readline()  # skip the first 3 lines, which includes the openbabel line, which changes every time
            output_str = file.read()

        expected_openbabel24 = '''7  6  0  0  0  0  0  0  0  0999 V2000
   -1.6424   -0.8401    0.0596 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4643   -1.2605    0.0894 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3236    0.1655   -0.0117 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3452   -0.1766    0.0125 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3469    0.1774   -0.0126 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5990    0.8179   -0.0580 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1824    1.1163   -0.0792 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  4  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  6  7  2  0  0  0  0
M  END

$$$$
'''
        expected_openbabel30 = '''7  6  0  0  0  0  0  0  0  0999 V2000
   -1.6424   -0.8401    0.0596 C   0  0  0  0  0  2  0  0  0  0  0  0
   -2.4643   -1.2605    0.0894 C   0  0  0  0  0  1  0  0  0  0  0  0
    0.3236    0.1655   -0.0117 C   0  0  0  0  0  1  0  0  0  0  0  0
   -0.3452   -0.1766    0.0125 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3469    0.1774   -0.0126 C   0  0  0  0  0  3  0  0  0  0  0  0
    1.5990    0.8179   -0.0580 C   0  0  0  0  0  3  0  0  0  0  0  0
    2.1824    1.1163   -0.0792 C   0  0  0  0  0  2  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  4  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  6  7  2  0  0  0  0
M  END

$$$$
'''

        assert output_str.strip() == expected_openbabel24.strip() or output_str.strip() == expected_openbabel30.strip()

    def test_select_atoms_missing_atoms(self):
        cmd = "exact c2 --input 4-helicene.mol --select-atoms 15-19,1,2,50"
        self.run_args(cmd)

        with open(os.path.join(self.results_folder, "extra.txt"), 'r') as file:
            output_str = file.read()
        
        message = "SELECT-ATOMS: 8 values selceted but only 7 of these values exist in molecule. Calculating using existing values."
        assert message in output_str

        cmd = "exact c2 --input 4-helicene.mol --select-atoms 50-60"
        dictionary_args = get_parsed_args(cmd.split())
        dictionary_args["argument_string"] = cmd + "\n"
        try:
            calc(dictionary_args) # run calc instead csm_run, because csm_run exits when calc raises exception
        except ValueError as err:
            assert str(err) == "select-atoms values do not exist in molecule."
        else:
            assert False

    def test_select_chains(self):
        cmd = "approx c2 --input 4yu4-protein.pdb --use-sequence --use-backbone --use-chains --select-chains A,C"
        results = self.run_args(cmd)
        assert results[0][0].csm == pytest.approx(0.006558605239814774, rel=1e-8)

        cmd = "approx c2 --input 4yu4-protein.pdb --use-sequence --use-backbone --use-chains --select-chains B,D"
        results = self.run_args(cmd)
        assert results[0][0].csm == pytest.approx(0.00821115474720635, rel=1e-8)

        cmd = "approx c2 --input 4yu4-protein.pdb --use-sequence --use-chains --select-chains A,C"
        results = self.run_args(cmd)
        assert results[0][0].csm == pytest.approx(0.020046822531283315, rel=1e-8)

        cmd = "approx c2 --input 4yu4-protein.pdb --use-sequence --use-chains --select-chains B,D"
        results = self.run_args(cmd)
        assert results[0][0].csm == pytest.approx(0.030496281388192603, rel=1e-8)

    def test_select_chains_missing_chains(self):
        cmd = "approx c2 --input 4yu4-protein.pdb --use-sequence --use-chains --select-chains B,E"
        dictionary_args = get_parsed_args(cmd.split())
        dictionary_args["argument_string"] = cmd + "\n"
        try:
            calc(dictionary_args) # run calc instead csm_run, because csm_run exits when calc raises exception
        except ValueError as err:
            assert str(err) == "select-chains values do not exist in molecule."
        else:
            assert False

        cmd = "approx c2 --input 4yu4-protein.pdb --use-sequence --use-chains --select-chains F-H"
        dictionary_args = get_parsed_args(cmd.split())
        dictionary_args["argument_string"] = cmd + "\n"
        try:
            calc(dictionary_args) # run calc instead csm_run, because csm_run exits when calc raises exception
        except ValueError as err:
            assert str(err) == "select-chains values do not exist in molecule."
        else:
            assert False

    def test_select_chains_on_dir(self):
        cmd = "approx c2 --input pdb-dir_BC --use-sequence --use-chains"
        reference = self.run_args(cmd)

        cmd = "approx c2 --input pdb-dir --use-sequence --use-chains --select-chains B,C"
        results = self.run_args(cmd)

        assert (len(results) == len(reference))

        for res, ref in zip(results, reference):
            assert res[0].csm == ref[0].csm

    def test_select_res(self):
        cmd0 = "approx c3 --input 7to4.pdb --use-sequence --select-res 15-306"
        results0 = self.run_args(cmd0)
        assert results0[0][0].csm == pytest.approx(0.11870267111205868, rel=1e-8)

        cmd1 = "approx c3 --input 7to4.pdb --use-sequence --select-res 330-530"
        results1 = self.run_args(cmd1)
        assert results1[0][0].csm == pytest.approx(12.508577670669041, rel=1e-8)

        cmd2 = "approx c3 --input 7to4.pdb --use-sequence --select-res 15-626"
        results2 = self.run_args(cmd2)
        assert results2[0][0].csm == pytest.approx(2.825056239763979, rel=1e-8)

        cmd3 = "approx c3 --input 7to4.pdb --use-sequence --select-res 687-1148"
        results3 = self.run_args(cmd3)
        assert results3[0][0].csm == pytest.approx(0.02553889330890735, rel=1e-8)

    def test_select_res_missing_res(self):
        cmd = "approx c3 --input 7to4.pdb --use-sequence --select-res 687-2000"
        self.run_args(cmd)

        with open(os.path.join(self.results_folder, "extra.txt"), 'r') as file:
            output_str = file.read()
        
        message = "SELECT-RES: 1314 values selceted but only 476 of these values exist in molecule. Calculating using existing values."
        assert message in output_str

        cmd = "approx c3 --input 7to4.pdb --use-sequence --select-res 2000-2500"
        dictionary_args = get_parsed_args(cmd.split())
        dictionary_args["argument_string"] = cmd + "\n"
        try:
            calc(dictionary_args) # run calc instead csm_run, because csm_run exits when calc raises exception
        except ValueError as err:
            assert str(err) == "select-res values do not exist in molecule."
        else:
            assert False
        

    def test_select_atoms_remove_hy(self):
        # --select-atoms removes specific atoms.
        # --remove-hy removes 'H' atoms.

        cmd = "exact c2 --input 4-helicene.mol --select-atoms 15-20,1-3 --remove-hy"
        results = self.run_args(cmd)
        assert len(results[0][0].molecule) == 7

        # test output
        with open(os.path.join(self.results_folder, "resulting_symmetric_coordinates.mol"), 'r') as file:
            file.readline()
            file.readline()
            file.readline()  # skip the first 3 lines, which includes the openbabel line, which changes every time
            output_str = file.read()

        expected_openbabel24 = '''7  6  0  0  0  0  0  0  0  0999 V2000
   -1.6424   -0.8401    0.0596 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4643   -1.2605    0.0894 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3236    0.1655   -0.0117 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3452   -0.1766    0.0125 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3469    0.1774   -0.0126 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5990    0.8179   -0.0580 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1824    1.1163   -0.0792 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  4  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  6  7  2  0  0  0  0
M  END

$$$$
'''
        expected_openbabel30 = '''7  6  0  0  0  0  0  0  0  0999 V2000
   -1.6424   -0.8401    0.0596 C   0  0  0  0  0  2  0  0  0  0  0  0
   -2.4643   -1.2605    0.0894 C   0  0  0  0  0  1  0  0  0  0  0  0
    0.3236    0.1655   -0.0117 C   0  0  0  0  0  1  0  0  0  0  0  0
   -0.3452   -0.1766    0.0125 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3469    0.1774   -0.0126 C   0  0  0  0  0  3  0  0  0  0  0  0
    1.5990    0.8179   -0.0580 C   0  0  0  0  0  3  0  0  0  0  0  0
    2.1824    1.1163   -0.0792 C   0  0  0  0  0  2  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  4  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  6  7  2  0  0  0  0
M  END

$$$$
'''

        assert output_str.strip() == expected_openbabel24.strip() or output_str.strip() == expected_openbabel30.strip()

    def test_use_backbone(self):
        cmd = "approx c2 --use-backbone --input 3alb-gkt4-h.pdb "
        result = self.run_args(cmd)
        assert len(result[0][0].molecule) == 48

        cmd = "approx c3 --input 2RLA-s3.pdb --use-backbone"
        result = self.run_args(cmd)
        assert len(result[0][0].molecule) == 12

    def test_use_backbone_use_seq(self):
        cmd = "approx c2 --input 1m2d.pdb --use-sequence --use-chains --use-backbone"
        result = self.run_args(cmd)

        cmd = "approx c2 --input 1m2d_bb.pdb --use-sequence --use-chains --use-backbone"
        result_bb = self.run_args(cmd)

        assert len(result[0][0].molecule) == len(result_bb[0][0].molecule)
        assert result[0][0].csm == pytest.approx(result[0][0].csm, rel=1e-8)

        cmd = "comfile cmd_use_backbone.txt --input 1m2d.pdb"
        result = self.run_args(cmd)

        assert len(result[0][0].molecule) == len(result_bb[0][0].molecule)
        assert result[0][0].csm == pytest.approx(result[0][0].csm, rel=1e-8)

    def test_ignore_atoms(self):
        def strip(myString):
            myString = myString.replace(' ', '').replace('\t', '').replace('\n', '')
            return myString
        # --ignore-atoms removes specific atoms.
        # The test checks that the result is identical to the result by use --select-atoms

        # cmd = "exact c2 --input 4-helicene.mol --select-atoms 15-19,1,2"  # len(_atoms) = 30
        cmd = "exact c2 --input 4-helicene.mol --ignore-atoms 3-14,20-30"
        results = self.run_args(cmd)
        assert len(results[0][0].molecule) == 7

        # test output
        with open(os.path.join(self.results_folder, "resulting_symmetric_coordinates.mol"), 'r') as file:
            file.readline()
            file.readline()
            file.readline()  # skip the first 3 lines, which includes the openbabel line, which changes every time
            output_str = file.read()

            expected_openbabel24 = '''7  6  0  0  0  0  0  0  0  0999 V2000
           -1.6424   -0.8401    0.0596 C   0  0  0  0  0  0  0  0  0  0  0  0
           -2.4643   -1.2605    0.0894 C   0  0  0  0  0  0  0  0  0  0  0  0
            0.3236    0.1655   -0.0117 C   0  0  0  0  0  0  0  0  0  0  0  0
           -0.3452   -0.1766    0.0125 C   0  0  0  0  0  0  0  0  0  0  0  0
            0.3469    0.1774   -0.0126 C   0  0  0  0  0  0  0  0  0  0  0  0
            1.5990    0.8179   -0.0580 C   0  0  0  0  0  0  0  0  0  0  0  0
            2.1824    1.1163   -0.0792 C   0  0  0  0  0  0  0  0  0  0  0  0
          1  2  1  0  0  0  0
          1  4  1  0  0  0  0
          3  4  1  0  0  0  0
          4  5  2  0  0  0  0
          5  6  1  0  0  0  0
          6  7  2  0  0  0  0
        M  END

        $$$$
        '''
            expected_openbabel30 = '''7  6  0  0  0  0  0  0  0  0999 V2000
           -1.6424   -0.8401    0.0596 C   0  0  0  0  0  2  0  0  0  0  0  0
           -2.4643   -1.2605    0.0894 C   0  0  0  0  0  1  0  0  0  0  0  0
            0.3236    0.1655   -0.0117 C   0  0  0  0  0  1  0  0  0  0  0  0
           -0.3452   -0.1766    0.0125 C   0  0  0  0  0  0  0  0  0  0  0  0
            0.3469    0.1774   -0.0126 C   0  0  0  0  0  3  0  0  0  0  0  0
            1.5990    0.8179   -0.0580 C   0  0  0  0  0  3  0  0  0  0  0  0
            2.1824    1.1163   -0.0792 C   0  0  0  0  0  2  0  0  0  0  0  0
          1  2  1  0  0  0  0
          1  4  1  0  0  0  0
          3  4  1  0  0  0  0
          4  5  2  0  0  0  0
          5  6  1  0  0  0  0
          6  7  2  0  0  0  0
        M  END

        $$$$
        '''


        assert strip(output_str) == strip(expected_openbabel24) or strip(output_str) == strip(expected_openbabel30)

    def test_ignore_sym(self):
        # --ignore-sym ignores atomic symbols
        cmd = "comfile ignore-sym/testcommand.txt --input ignore-sym"
        results = self.run_args(cmd)
        assert len(results[0][0].molecule.equivalence_classes) == 2
        assert len(results[1][0].molecule.equivalence_classes) == 1
        assert len(results[2][0].molecule.equivalence_classes) == 2
        assert len(results[3][0].molecule.equivalence_classes) == 1

        # TODO add output test?

    def test_use_mass(self):
        # --use-mass uses atomic masses to define center of mass
        cmd = "trivial s8 --input ferrocene.xyz"
        result1 = self.run_args(cmd)
        cmd = "trivial s8 --input ferrocene.xyz --use-mass"
        results = self.run_args(cmd)
        assert result1[0][0].molecule.center_of_mass != results[0][0].molecule.center_of_mass

    def test_babel_bond(self):
        # --babel-bond computes bonding
        cmd = "trivial ch --input ferrocene.xyz"
        result1 = self.run_args(cmd)
        cmd = "trivial ch --input ferrocene.xyz --babel-bond"
        results = self.run_args(cmd)
        assert len(result1[0][0].molecule.bondset) != len(results[0][0].molecule.bondset)

    def test_use_sequence(self):
        # --use-sequence creates equivalence class based on sequence
        cmd = "approx c2 --use-sequence --input lig-4kem.pdb"
        result = self.run_args(cmd)
        assert result[0][0].csm == pytest.approx(0.0062574, abs=1e-5)

    def test_select_mols(self):
        # --select-mols when given a folder or multi-mol file
        # folder:
        cmd = "trivial c2 --input ignore-sym"
        results = self.run_args(cmd)
        assert len(results) == 4

        cmd = "trivial c2 --input ignore-sym --select-mols 1,3-4"
        results = self.run_args(cmd)
        assert len(results) == 3

        # file:
        cmd = "trivial c2 --input squarate.xyz"
        results = self.run_args(cmd)
        assert len(results) == 5

        cmd = "trivial c2 --input squarate.xyz --select-mols 1,3-4"
        results = self.run_args(cmd)
        assert len(results) == 3

        #same deals but with comfile
        cmd = "comfile comfile.txt --input squarate.xyz --select-mols 1,3-4"
        results = self.run_args(cmd)
        assert len(results) == 3

        cmd = "comfile comfile.txt --input ignore-sym --select-mols 1,3-4"
        results = self.run_args(cmd)
        assert len(results) == 3


    # output
    def test_legacy(self):
        # why is it printing filename instead of index all of a sudden?
        # --legacy-output prints old style csm format
        cmd = "exact c2 --input squarate.xyz --select-mols 1 --legacy-output --output {}/legacy.txt".format(
            self.results_folder)
        csm_run(cmd.split())
        with open(os.path.join(self.results_folder, "legacy.txt"), 'r') as ofile:
            out = ofile.read()
        with open(os.path.join("expected", "legacy.txt"), 'r') as ofile:
            exp = ofile.read()
        assert out.strip() == exp.strip()

    def test_json_output(self):
        # same problem as legacy
        # --json-output. only works with --legacy-output
        cmd = "exact c2 --input squarate.xyz  --json-output --select-mols 1,2"
        self.run_args(cmd)
        with open(os.path.join(self.results_folder, "json-results.json"), 'r') as ofile:
            out = json.loads(ofile.read())

    def test_legacy_files(self):
        cmd = "exact c2 --input squarate.xyz --legacy-files"
        self.run_args(cmd)
        # --legacy-files toggles whether we create a folder of legacy files or not
        output_path = os.path.join(self.results_folder, "old-csm-output")
        assert os.path.isdir(output_path)
        with open(os.path.join(output_path, "0001_L01_c2.xyz"), "r") as file:
            out = file.read()
        with open(os.path.join("expected", "example-legacy-file.xyz"), "r") as file:
            exp = file.read()
        assert out.strip() == exp.strip()

    def test_verbose(self):
        # verbose-- approx only
        cmd = "approx c2 --input squarate.xyz --verbose"
        self.run_args(cmd)
        output_path = os.path.join(self.results_folder, "approx")
        assert os.path.isdir(output_path)
        # todo: because the output contains a variable runtime, running a comparison is a bit tedious, leaving it for now


    # parallel
    def test_parallel_mols_in_file(self):
        cmd = "comfile comfile.txt --input many-mols.xyz --parallel"
        result1 = self.run_args(cmd)
        with open(os.path.join(self.results_folder, "csm.txt"), "r") as file:
            out = file.read()
        with open(os.path.join("parallel-out", "csm.txt"), "r") as file:
            exp = file.read()
        assert out == exp

    def test_sn_max(self):
        # --sn-max (relevant only for chirality)
        cmd = "exact ch --input bis(dth)copper(I).mol"
        result1 = self.run_args(cmd)
        cmd = "exact ch --input bis(dth)copper(I).mol --sn-max 2"
        result2 = self.run_args(cmd)
        assert result1[0][0].op_type == 'SN' and result1[0][0].op_order == 4
        assert result2[0][0].op_type == 'CS' and result2[0][0].op_order == 2

    def test_normalize(self):
        test_file = os.path.join('reading-fragments', 'water-6.mol')
        cmd = "approx c3 --fibonacci 1 --read-fragments --input " + test_file + " --normalize 0 1 2 3 4 5 6"
        results = self.run_args(cmd)
        hi = 1

    def test_format(self):
        # --in-format
        cmd = "exact cs --input 4-helicene.ttxt --keep-structure --in-format mol"
        results = self.run_args(cmd)
        assert results[0][0].csm == pytest.approx(0.793551, abs=1e-5)
        # --out-format
        cmd = "exact cs --input 4-helicene.mol --keep-structure --out-format bs"
        results = self.run_args(cmd)
        with open(os.path.join(self.results_folder, "resulting_symmetric_coordinates.bs"), 'r') as ofile:
            out = ofile.read()
        with open(os.path.join("expected", "testformatout.bs"), 'r') as efile:
            exp = efile.read()
        assert out == exp

    # comfile

    def test_old_cmd(self):
        # old-cmd in comfile
        cmd = "comfile comfile.txt --input reading-fragments"
        results1 = self.run_args(cmd)
        cmd = "comfile oldcomfile.txt --input reading-fragments --old-cmd"
        results2 = self.run_args(cmd)
        for r1, r2 in zip(results1, results2):
            for rr1, rr2 in zip(r1, r2):
                print(rr1.csm, rr2.csm)
                print(round(rr1.csm), round(rr2.csm))
                assert round(rr1.csm, 5) == round(rr2.csm, 5)

    # exact
    def test_use_perm(self):
        cmd = "exact c4 --input squarate.xyz --use-perm squarateperm.txt"
        results = self.run_args(cmd)
        assert results[0][0].perm == [0, 1, 2, 3, 4, 5, 6, 7]
        cmd = "exact cs --input cryptands-no-metal.sdf --use-perm perm_select_atoms.txt" \
              " --select-atoms 7,11-14 --select-mol 1"
        results = self.run_args(cmd)
        assert results[0][0].perm == [4, 3, 2, 1, 0]

    def test_keep_structure(self):
        cmd = "exact cs --input 4-helicene.mol --keep-structure"
        results = self.run_args(cmd)
        assert results[0][0].csm == pytest.approx(0.793551, abs=1e-5)

    def test_output_perms(self):
        try:
            with open(os.path.join(self.results_folder, "exact", "4-helicene_L01_cs.csv"), 'w') as file:
                # reset perms.csv
                pass
        except:
            pass
        cmd = "exact cs --input 4-helicene.mol --keep-structure --output-perms"
        self.run_args(cmd)
        out_rows = []
        with open(os.path.join(self.results_folder, "exact", "4-helicene_L01_cs.csv"), 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                out_rows.append(row)
        assert out_rows[0] == ['op', 'Permutation', 'Direction', 'CSM']
        assert out_rows[1][0]=='CS2'
        assert out_rows[2][1]=='[17, 29, 30, 27, 28, 25, 26, 23, 24, 18, 19, 20, 21, 22, 15, 16, 1, 10, 11, 12, 13, 14, 8, 9, 6, 7, 4, 5, 2, 3]'

    # approx
    def test_parallel_dirs(self):
        cmd = "approx c2 --input 3alb-gkt4-h.pdb --use-sequence --parallel-dirs"
        results = self.run_args(cmd)
        assert not results[0][0].failed

    def test_input_chain_perm(self):
        # cmd = "approx c2 --input 3alb-gkt4-h.pdb --input-chain-perm --verbose"
        # results = self.run_args(cmd)
        # assert results[0][0].chain_perm == [2, 1, 0, 3]  # when run with use-chains it gets 0 1 2 3

        cmd = "approx c3 --input 1nc7.pdb  --use-sequence --use-chains --input-chain-perm CBDA2130.txt"
        results = self.run_args(cmd)
        print("result:", results[0][0].chain_perm)
        assert results[0][0].chain_perm == [2,1,3,0]

        # cmd = "approx c4 --input 1v0z.pdb --use-sequence --input-chain-perm BCDA1230.txt"
        # results = self.run_args(cmd)
        # print("result:", results[0][0].chain_perm)
        # assert results[0][0].chain_perm == [1,2,3,0]

        # cmd = "approx c4 --input 1v0z.pdb --use-sequence --input-chain-perm DABC3012.txt"
        # results = self.run_args(cmd)
        # print("result:", results[0][0].chain_perm)
        # assert results[0][0].chain_perm == [3,0,1,2]

    def test_input_chain_perm_trivial(self):
        cmd = "trivial c2 --input 3alb-gkt4-h.pdb --input-chain-perm --verbose"
        results = self.run_args(cmd)
        assert results[0][0].chain_perm == [2, 1, 0, 3]
        with open(os.path.join(self.results_folder, "trivial", "3alb-gkt4-h_L01_c2.tsv"), 'r') as tsvfile:
            pass

    def test_unequal_chains(self):
        cmd = "approx c2 --input 4a5k-467.pdb --use-chains"
        results = self.run_args(cmd)
        assert results[0][0].chain_perm==[0,3,2,1]

        cmd = "trivial c2 --input 4a5k-467.pdb --use-chains"
        results = self.run_args(cmd)
        assert results[0][0].chain_perm==[0,3,2,1]
    

    def test_use_chains(self):
        cmd = "approx c3 --input 2RLA-s3.pdb --use-chains"
        results = self.run_args(cmd)
        assert len(results[0][0].molecule.chains) == 3
        cmd = "approx c3 --input 2RLA-s3.pdb"
        results = self.run_args(cmd)
        assert len(results[0][0].molecule.chains) == 1

    def test_selective(self):
        cmd = "approx c3 --input 2RLA-s3.pdb --fibonacci 50 --selective 3"
        results = self.run_args(cmd)

    def test_permute_chains(self):
        cmd = "trivial c3 --input 2RLA-s3.pdb"
        results1 = self.run_args(cmd)
        assert results1[0][0].csm == pytest.approx(49.078986, abs=1e-5)

        cmd = "trivial c3 --input 2RLA-s3.pdb --permute-chains"
        results2 = self.run_args(cmd)
        assert results2[0][0].csm == pytest.approx(0.009663, abs=1e-5)

    def test_directions(self):
        # TODO: refine detect-outliers test
        # default
        cmd = "approx c3 --input 3alb-gkt4-h.pdb"
        results = self.run_args(cmd)
        assert len(results[0][0].ongoing_statistics['approx']['c3']) == 9
        # detect-outliers
        cmd = "approx c3 --input 3alb-gkt4-h.pdb --detect-outliers"
        results = self.run_args(cmd)
        assert len(results[0][0].ongoing_statistics['approx']['c3']) == 18
        # no-orthogonal
        cmd = "approx c3 --input 2RLA-s3.pdb --no-orthogonal"
        results = self.run_args(cmd)
        assert len(results[0][0].ongoing_statistics['approx']['c3']) == 3
        # use-best-dir
        cmd = "approx c3 --input 2RLA-s3.pdb --use-best-dir"
        results = self.run_args(cmd)
        assert len(results[0][0].ongoing_statistics['approx']['c3']) == 3
        # use-best-dir + --no-orthogonal
        cmd = "approx c3 --input 2RLA-s3.pdb --use-best-dir --no-orthogonal"
        results = self.run_args(cmd)
        assert len(results[0][0].ongoing_statistics['approx']['c3']) == 1
        # fibonacci
        cmd = "approx c3 --input 2RLA-s3.pdb --fibonacci 20"
        results = self.run_args(cmd)
        assert len(results[0][0].ongoing_statistics['approx']['c3']) == 20
        # dir
        cmd = "approx c3 --input 2RLA-s3.pdb --dir 1 0 0"
        results = self.run_args(cmd)
        assert len(results[0][0].ongoing_statistics['approx']['c3']) == 1

    def test_algorithms(self):
        # default
        cmd = "approx c3 --input 2RLA-s3.pdb"
        results = self.run_args(cmd)
        assert round(results[0][0].csm, 6) == 0.009663
        stats1 = tuple(results[0][0].ongoing_statistics['approx']['c3'][3]['stats']['dirs'][0])
        # greedy
        cmd = "approx c3 --input 2RLA-s3.pdb --greedy"
        results = self.run_args(cmd)
        assert round(results[0][0].csm, 6) == 0.009663
        stats2 = tuple(results[0][0].ongoing_statistics['approx']['c3'][3]['stats']['dirs'][0])
        # many-chains
        cmd = "approx c3 --input 2RLA-s3.pdb --many-chains"
        results = self.run_args(cmd)
        assert round(results[0][0].csm, 6) == 0.009663
        stats3 = tuple(results[0][0].ongoing_statistics['approx']['c3'][3]['stats']['dirs'][0])
        # keep-structure
        cmd = "approx c3 --input 2RLA-s3.pdb --keep-structure"
        results = self.run_args(cmd)
        assert round(results[0][0].csm, 6) == 0.009663
        stats4 = tuple(results[0][0].ongoing_statistics['approx']['c3'][3]['stats']['dirs'][0])

        test = set([stats1, stats2, stats3, stats4])
        assert len(test) == 4


class TestFragments(RunThings):
    os.chdir(test_dir)
    results_folder = "csm_tests"

    # shutil.rmtree(results_folder)
    # os.mkdir(results_folder)

    def run_args(self, args_str):
        return super()._run_args(args_str, self.results_folder)

    def test_read_fragments_baseline(self):
        # --read-fragments reads from .mol, .pdb as chains
        filenames = ["model-endmdl-withIDS.pdb",
                     #"model-endmdl-withoutIDS.pdb",
                     "water-6.mol",]
                     #"water-6.pdb"]

        for filename in filenames:
            cmd = "approx c3 --fibonacci 1 --input "
            dir_file_name = os.path.join("reading-fragments", filename)
            cmd += dir_file_name
            results = self.run_args(cmd)
            print(cmd, len(results[0][0].molecule.chains))
            assert len(results[0][0].molecule.chains) == 1

    def test_read_fragments_with_read_fragments(self):
        # --read-fragments reads from .mol, .pdb as chains
        filenames = [
            "model-endmdl-withIDS.pdb",
            "model-endmdl-withoutIDS.pdb",
            "water-6.mol",
            "water-6.pdb"]
        for filename in filenames:
            cmd = "approx c3 --fibonacci 1 --read-fragments --input "
            dir_file_name = os.path.join("reading-fragments", filename)
            cmd += dir_file_name
            results = self.run_args(cmd)
            print(cmd)
            assert len(results[0][0].molecule.chains) == 6

    def test_read_fragments_with_read_fragments_and_use_chain(self):
        # --read-fragments reads from .mol, .pdb as chains
        filenames = ["model-endmdl-withIDS.pdb",
                     "model-endmdl-withoutIDS.pdb",
                     "water-6.mol",
                     "water-6.pdb"]
        for filename in filenames:
            cmd = "approx c3 --fibonacci 1 --read-fragments --use-chains --input "
            dir_file_name = os.path.join("reading-fragments", filename)
            cmd += dir_file_name
            results = self.run_args(cmd)
            print(cmd)
            assert len(results[0][0].molecule.chains) == 6


class TestChirality(RunThings):
    results_folder = "csm_tests"

    def run_args(self, args_str):
        return super()._run_args(args_str, self.results_folder)

    def test_trivial(self):
        # --babel-bond computes bonding
        cmd = "trivial ch --input ferrocene.xyz"
        result1 = self.run_args(cmd)
        assert result1[0][0].csm == pytest.approx(17.152493, abs=1e-5)
        assert result1[0][0].overall_statistics["best chirality"] == "CS"

    def test_trivial_2(self):
        cmd = "trivial ch --input ferrocene.xyz --babel-bond"
        results = self.run_args(cmd)
        assert results[0][0].csm == pytest.approx(17.152493, abs=1e-5)
        assert results[0][0].overall_statistics["best chirality"] == "CS"

    def test_trivial_3(self):
        cmd = "trivial ch --input 2RLA-s3.pdb"
        result1 = self.run_args(cmd)
        assert result1[0][0].csm == pytest.approx(0.111737, abs=1e-5)
        assert result1[0][0].overall_statistics["best chirality"] == "CS"

    def test_trivial_4(self):
        cmd = "trivial ch --input 2RLA-s3.pdb --permute-chains"
        result2 = self.run_args(cmd)
        assert result2[0][0].csm == pytest.approx(0.111737, abs=1e-5)
        assert result2[0][0].overall_statistics["best chirality"] == "CS"

    def test_exact(self):
        # --sn-max (relevant only for chirality)
        cmd = "exact ch --input bis(dth)copper(I).mol"  # matches
        result1 = self.run_args(cmd)
        assert result1[0][0].csm == pytest.approx(0, abs=1e-5)
        assert result1[0][0].overall_statistics["best chirality"] == "S4"

    def test_exact_2(self):
        cmd = "exact ch --input bis(dth)copper(I).mol --sn-max 2"  # matches
        result2 = self.run_args(cmd)
        assert result2[0][0].csm == pytest.approx(4.394381, abs=1e-5)
        assert result2[0][0].overall_statistics["best chirality"] == "CS"

    def test_approx(self):
        cmd = "approx ch --input ferrocene.xyz --babel-bond"  # matches
        results = self.run_args(cmd)
        assert results[0][0].csm == pytest.approx(0.241361, abs=1e-5)
        assert results[0][0].overall_statistics["best chirality"] == "CS"

    def test_approx_2(self):
        cmd = "approx ch --input bis(dth)copper(I).mol"  # matches
        result1 = self.run_args(cmd)
        assert result1[0][0].csm == pytest.approx(0, abs=1e-5)
        assert result1[0][0].overall_statistics["best chirality"] == "S4"

    def test_approx_3(self):
        cmd = "approx ch --input bis(dth)copper(I).mol --sn-max 2"  # matches
        result2 = self.run_args(cmd)
        assert result2[0][0].csm == pytest.approx(4.394381, abs=1e-5)
        assert result2[0][0].overall_statistics["best chirality"] == "CS"


class xTestComplicated(RunThings):
    def __init__(self):
        # an option that is too complicated to test and will simply have to be checked by eye:
        # 1. the default of writing results to a file with a timestamp (rather than with overwrite)
        # 2. --pipe
        # options that are difficult to test that they're working as expected:
        # 3. --parallel
        # 4. --parallel-approx
        # 5. print-approx
        # 6. simple
        # 7. timeout, global-timeout
        pass
