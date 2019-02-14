'''
A set of tests to check that a variety of input arguments are still valid and haven't broken from changes introduced to the code-- doesn't test validity of result, just
that the code runs successfully and returns A result.
'''
import csv

import os
import pytest
import shutil

from csm.main.csm_run import csm_run


class Runner:
    def run_args(self, args_str):
        args = args_str.split()
        results_arr = csm_run(args)
        result = results_arr[0][0]
        assert result.csm is not None


class xTestExact(Runner):
    def xtest_plain(self):
        self.run_args(
            r"exact c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt")

    @pytest.mark.parametrize("run_str", [
        r"exact ch --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --sn-max 6",
        r"exact c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2RLA-3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --remove-hy",
        r"exact c3 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --use-mass",
        r"exact s8 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --ignore-sym",
        r"exact s6 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --read-fragments",
        r"exact s2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2RLA-3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --use-sequence",
        r"exact ci --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --use-chains",
        r"exact c4 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --babel-bond",
    ])
    def xtest_input_args(self, run_str):
        self.run_args(run_str)

    def xtest_exact_args(self):
        self.run_args(
            r"exact c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --keep-structure")

    @pytest.mark.parametrize("run_str",
                             [
                                 r"exact c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --print-denorm --legacy",
                                 r"exact c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --output-branches --legacy",
                                 r"exact c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --output-perms perms.csv --legacy",
                                 r"exact c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --print-local --legacy",
                                 r"exact c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --json-output --legacy",
                                 r"exact c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt --simple",
                             ])
    def xtest_output_args(self, run_str):
        self.run_args(run_str)


class xTestApprox(Runner):
    def xtest_plain(self):
        self.run_args(
            r"approx c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt ")
        self.run_args(
            r"approx c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt  --print-approx")

    @pytest.mark.parametrize("run_str",
                             [

                                 r"approx ch --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt  --sn-max 6",
                                 r"approx ci --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt  --use-chains",
                                 r"approx s2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt  --use-sequence",
                                 r"approx s6 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt  --read-fragments",

                             ])
    def xtest_input_args(self, run_str):
        self.run_args(run_str)

    @pytest.mark.parametrize("run_str",
                             [
                                 r"approx c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\approx\csm_testoutput.txt  --detect-outliers",
                                 r"approx c3 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\approx\csm_testoutput.txt  --no-orthogonal",
                                 r"approx c5 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\approx\csm_testoutput.txt  --detect-outliers --no-orthogonal",
                                 r"approx cs --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\approx\csm_testoutput.txt  --fibonacci 20 --parallel 2",
                                 r"approx c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\approx\csm_testoutput.txt  --use-best-dir",
                                 r"approx c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\approx\csm_testoutput.txt  --many-chains",
                                 r"approx c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\approx\csm_testoutput.txt  --greedy",
                                 r"approx c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\approx\csm_testoutput.txt  --fibonacci 20 --selective 3",
                                 r"approx c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\approx\csm_testoutput.txt  --polar",
                             ])
    def xtest_approx_args(self, run_str):
        self.run_args(run_str)

    @pytest.mark.parametrize("run_str",
                             [
                                 r"approx c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\approx --output C:\Users\devora\Sources\temp\csm_argument_tests\approx\multi  --use-best-dir",
                                 r"approx c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\approx --output C:\Users\devora\Sources\temp\csm_argument_tests\approx\multi  --fibonacci 3",
                             ])
    def xtest_multi_mols(self, run_str):
        self.run_args(run_str)

    def xtest_output_args(self):
        self.run_args(
            r"approx c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt  --json-output --legacy")


class xTestTrivial(Runner):
    def xtest_plain(self):
        self.run_args(
            r"trivial c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\just-one-mol.mol --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt ")

    @pytest.mark.parametrize("run_str", [

        r"trivial ch --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt  --sn-max 6",
        r"trivial c4 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt  --babel-bond",
        r"trivial ci --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt  --use-chains",
        r"trivial s2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt  --use-sequence",
        r"trivial s6 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests\2m7w-q3.pdb --output C:\Users\devora\Sources\temp\csm_argument_tests\csm_testoutput.txt  --read-fragments",
    ])
    def xtest_input_args(self, run_str):
        self.run_args(run_str)

    @pytest.mark.parametrize("run_str", [

        r"trivial ch --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests --output C:\Users\devora\Sources\temp\csm_argument_tests\multi_molecules  --sn-max 6",
        r"trivial c4 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests --output C:\Users\devora\Sources\temp\csm_argument_tests\multi_molecules  --babel-bond",
        r"trivial ci --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests --output C:\Users\devora\Sources\temp\csm_argument_tests\multi_molecules  --use-chains",
        r"trivial s2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests --output C:\Users\devora\Sources\temp\csm_argument_tests\multi_molecules  --use-sequence",
        r"trivial s6 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests --output C:\Users\devora\Sources\temp\csm_argument_tests\multi_molecules  --read-fragments",
        r"trivial c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests --output C:\Users\devora\Sources\temp\csm_argument_tests\multi_molecules  --select-atoms 1,3,4-6",
        r"trivial c2 --input C:\Users\devora\Sources\csm\python\tests\unit_tests\molecules_for_tests --output C:\Users\devora\Sources\temp\csm_argument_tests\multi_molecules  --select-mols 1,2",

    ])
    def xtest_multiple_mols(self, run_str):
        self.run_args(run_str)


class RunThings():
    def run_args(self, args_str):
        args_str += " --output csm_tests --overwrite"
        args = args_str.split()
        results_arr = csm_run(args)
        return results_arr


class TestBasic(RunThings):
    test_dir = r"C:\Users\devora\Sources\csm\csm\python\tests\argument_tests\files_for_tests"
    os.chdir(test_dir)
    shutil.rmtree("csm_tests")
    os.mkdir("csm_tests")

    # an option that is too complicated to test and will simply have to be checked by eye:
    # 1. the default of writing results to a file with a timestamp (rather than with overwrite)
    # 2. --pipe
    # options that are difficult to test that they're working as expected:
    # 3. --parallel
    # 4. --parallel-approx
    # 5. print-approx
    # 6. simple
    # 7. timeout, global-timeout

    # input
    def xtest_connect(self):
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

    def xtest_remove_hy(self):
        # --remove-hy removes hy. can't be used with select-atoms

        # baseline:
        # cmd="exact cs --input 4-helicene.mol --keep-structure"
        # results=self.run_args(cmd)
        # assert len(results[0][0].molecule) == 30

        cmd = "exact cs --input 4-helicene.mol --remove-hy"
        results = self.run_args(cmd)
        assert len(results[0][0].molecule) == 18

        # test output
        with open("csm_tests\\resulting_symmetric_coordinates.mol", 'r') as file:
            file.readline()
            file.readline()  # skip the openbabel line, which changes every time
            output_str = file.read()
        with open("removehyexpected.mol", 'r') as file:
            file.readline()
            file.readline()  # skip the openbabel line, which changes every time
            expected_str = file.read()

        assert expected_str == output_str

    def xtest_select_atoms(self):
        # --select-atoms removes sepcific atoms. cannot be used in conjunction with remove-hy
        # baseline:
        # cmd="exact c2 --input 4-helicene.mol --keep-structure"
        # results=self.run_args(cmd)
        # assert len(results[0][0].molecule) == 30

        cmd = "exact c2 --input 4-helicene.mol --select-atoms 15-19,1,2"
        results = self.run_args(cmd)
        assert len(results[0][0].molecule) == 7

        # test output
        with open("csm_tests\\resulting_symmetric_coordinates.mol", 'r') as file:
            file.readline()
            file.readline()  # skip the openbabel line, which changes every time
            output_str = file.read()
        with open("selectatomsexpected.mol", 'r') as file:
            file.readline()
            file.readline()  # skip the openbabel line, which changes every time
            expected_str = file.read()

        assert expected_str == output_str

    def xtest_ignore_sym(self):
        # --ignore-sym ignores atomic symbols
        cmd = "comfile ignore-sym/testcommand.txt --input ignore-sym"
        results = self.run_args(cmd)
        assert len(results[0][0].molecule.equivalence_classes) == 2
        assert len(results[1][0].molecule.equivalence_classes) == 1
        assert len(results[2][0].molecule.equivalence_classes) == 2
        assert len(results[3][0].molecule.equivalence_classes) == 1

        # TODO add output test?

    def xtest_use_mass(self):
        # --use-mass uses atomic masses to define center of mass
        cmd = "trivial s8 --input ferrocene.xyz"
        result1 = self.run_args(cmd)
        cmd = "trivial s8 --input ferrocene.xyz --use-mass"
        results = self.run_args(cmd)
        assert result1[0][0].molecule.center_of_mass != results[0][0].molecule.center_of_mass

    def xtest_babel_bond(self):
        # --babel-bond computes bonding
        cmd = "trivial ch --input ferrocene.xyz"
        result1 = self.run_args(cmd)
        cmd = "trivial ch --input ferrocene.xyz --babel-bond"
        results = self.run_args(cmd)
        assert len(result1[0][0].molecule.bondset) != len(results[0][0].molecule.bondset)

    def xtest_use_sequence(self):
        # --use-sequence creates equivalence class based on sequence
        cmd = "approx c2 --use-sequence --input lig-4kem.pdb"
        result = self.run_args(cmd)
        assert result[0][0].csm == pytest.approx(0.0062574, abs=1e-5)

    def todo_test_read_fragments(self):
        # --read-fragments reads from .mol, .pdb as chains
        assert False

    def xtest_select_mols(self):
        # --select-mols when given a folder or multi-mol file
        #folder:
        cmd = "trivial c2 --input ignore-sym"
        results = self.run_args(cmd)
        assert len(results) ==4

        cmd = "trivial c2 --input ignore-sym --select-mols 1,3-4"
        results = self.run_args(cmd)
        assert len(results) ==3

        #file:
        cmd = "trivial c2 --input squarate.xyz"
        results = self.run_args(cmd)
        assert len(results) ==5

        cmd = "trivial c2 --input squarate.xyz --select-mols 1,3-4"
        results = self.run_args(cmd)
        assert len(results) ==3

        #TODO: output tests


    # output
    def xtest_legacy(self):
        # --legacy-output prints old style csm format
        cmd = "exact c2 --input squarate.xyz --select-mols 1 --legacy-output --output csm_tests/legacy.txt"
        csm_run(cmd.split())
        with open("csm_tests\legacy.txt", 'r') as ofile:
            out = ofile.read()
        with open("legacy.txt", 'r') as ofile:
            exp = ofile.read()
        assert out == exp

    def xtest_json_output(self):
        # --json-output. only works with --legacy-output
        cmd = "exact c2 --input squarate.xyz --select-mols 1 --legacy-output --output csm_tests/legacy.json --json-output"
        csm_run(cmd.split())
        with open("csm_tests\legacy.json", 'r') as ofile:
            out = ofile.read()
        with open("legacy.json", 'r') as ofile:
            exp = ofile.read()
        assert out == exp

    def xtest_print_denorm(self):
        # --print-denorm prints denormalized coordinates for original molecule instead of normalized
        cmd = "exact c2 --input squarate.xyz --select-mols 1 --print-denorm"
        self.run_args(cmd)
        with open("csm_tests/initial_normalized_coordinates.xyz", "r") as file:
            output2 = file.read()

        cmd = "exact c2 --input squarate.xyz --select-mols 1"
        self.run_args(cmd)
        with open("csm_tests/initial_normalized_coordinates.xyz", "r") as file:
            output1 = file.read()

        assert output1 != output2

    def todo_test_legacy_files(self):
        cmd = "exact c2 --input squarate.xyz --legacy-files"

        # --legacy-files toggles whether we create a folder of legacy files or not
        assert False

    def todo_test_verbose(self):
        # verbose-- approx only
        assert False

    # shared stuff

    def todo_test_sn_max(self):
        # --sn-max (relevant only for chirality)
        assert False

    def todo_test_normalize(self):
        # --normalize
        assert False

    def xtest_format(self):
        # --in-format
        cmd = "exact cs --input 4-helicene.ttxt --keep-structure --in-format mol"
        results = self.run_args(cmd)
        assert results[0][0].csm == pytest.approx(0.793551, abs=1e-5)
        # --out-format
        cmd = "exact cs --input 4-helicene.mol --keep-structure --out-format bs"
        results = self.run_args(cmd)
        with open("csm_tests/resulting_symmetric_coordinates.bs", 'r') as ofile:
            out = ofile.read()
        with open("testformatout.bs", 'r') as efile:
            exp = efile.read()
        assert out == exp

    def todo_test_old_cmd(self):
        # old-cmd in comfile
        assert False

    # exact
    def xtest_use_perm(self):
        cmd = "exact c4 --input squarate.xyz --use-perm squarateperm.txt"
        results = self.run_args(cmd)
        assert results[0][0].perm == [0, 1, 2, 3, 4, 5, 6, 7]

    def xtest_keep_structure(self):
        cmd = "exact cs --input 4-helicene.mol --keep-structure"
        results = self.run_args(cmd)
        assert results[0][0].csm == pytest.approx(0.793551, abs=1e-5)

    def xtest_output_perms(self):
        cmd = "exact cs --input 4-helicene.mol --keep-structure --output-perms csm_tests/perms.csv"
        self.run_args(cmd)
        out_rows = []
        with open("csm_tests/perms.csv", 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                out_rows.append(row)
        assert out_rows == [
            ['Permutation', 'Direction', 'CSM'],
            [
                '[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]',
                '[0.042475 0.020727 0.998883]', '0.7935529301987154'],
            [
                '[17, 29, 30, 27, 28, 25, 26, 23, 24, 18, 19, 20, 21, 22, 15, 16, 1, 10, 11, 12, 13, 14, 8, 9, 6, 7, 4, 5, 2, 3]',
                '[-0.514303 -0.856692  0.039646]', '0.7935519339113517']
        ]

    # approx
    def xtest_use_chains(self):
        cmd = "approx c3 --input 2rla-s3.pdb --use-chains"
        results = self.run_args(cmd)
        assert len(results[0][0].molecule.chains)==3
        cmd = "approx c3 --input 2rla-s3.pdb"
        results = self.run_args(cmd)
        assert len(results[0][0].molecule.chains)==1

    def todo_test_directions(self):
        # TODO: refine detect-outliers test
        # default
        cmd = "approx c3 --input crown21.xyz --connect crown21connectivity.txt"
        results = self.run_args(cmd)
        test1 = [item['dir'] for item in results[0][0].ongoing_statistics['approx']]
        # detect-outliers
        cmd = "approx c3 --input crown21.xyz --connect crown21connectivity.txt --detect-outliers"
        results = self.run_args(cmd)
        test2 = [item['dir'] for item in results[0][0].ongoing_statistics['approx']]
        # no-orthogonal
        cmd = "approx c3 --input 2rla-s3.pdb --no-orthogonal"
        results = self.run_args(cmd)
        assert len(results[0][0].ongoing_statistics['approx']) == 3
        # use-best-dir
        cmd = "approx c3 --input 2rla-s3.pdb --use-best-dir"
        results = self.run_args(cmd)
        assert len(results[0][0].ongoing_statistics['approx']) == 3
        # use-best-dir + --no-orthogonal
        cmd = "approx c3 --input 2rla-s3.pdb --use-best-dir --no-orthogonal"
        results = self.run_args(cmd)
        assert len(results[0][0].ongoing_statistics['approx']) == 1
        # fibonacci
        cmd = "approx c3 --input 2rla-s3.pdb --fibonacci 20"
        results = self.run_args(cmd)
        assert len(results[0][0].ongoing_statistics['approx']) == 20
        # dir
        cmd = "approx c3 --input 2rla-s3.pdb --dir 1 0 0"
        results = self.run_args(cmd)
        assert len(results[0][0].ongoing_statistics['approx']) == 1

    def todo_test_algorithms(self):
        # default
        cmd = "approx c3 --input 2rla-s3.pdb"
        results = self.run_args(cmd)
        print("!!!!!!!!!!", results[0][0].csm)
        # greedy
        cmd = "approx c3 --input 2rla-s3.pdb --greedy"
        results = self.run_args(cmd)
        print("!!!!!!!!!!", results[0][0].csm)
        # many-chains
        cmd = "approx c3 --input 4a5k-467.pdb --many-chains"
        results = self.run_args(cmd)
        print("!!!!!!!!!!", results[0][0].csm)
        # keep-structure
        cmd = "approx c3 --input 2rla-s3.pdb --keep-structure"
        results = self.run_args(cmd)
        print("!!!!!!!!!!", results[0][0].csm)

    def xtest_selective(self):
        cmd = "approx c3 --input 2rla-s3.pdb --fibonacci 50 --selective 3"
        results = self.run_args(cmd)

    def xtest_dont_permute_chains(self):
        cmd = "trivial c3 --input 2rla-s3.pdb"
        results1 = self.run_args(cmd)
        assert results1[0][0].csm == pytest.approx(0.009663, abs=1e-5)
        # dont-permute-chains-- deactivate use-chains
        cmd = "trivial c3 --input 2rla-s3.pdb --dont-permute-chains"
        results2 = self.run_args(cmd)
        assert results2[0][0].csm == pytest.approx(49.078986, abs=1e-5)
