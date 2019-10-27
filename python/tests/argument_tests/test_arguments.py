'''
A set of tests to check that a variety of input arguments are still valid and haven't broken from changes introduced to the code-- doesn't test validity of result, just
that the code runs successfully and returns A result.
'''
import csv

import os
import pytest
import shutil

from csm.main.csm_run import csm_run

class RunThings():
    def _run_args(self, args_str, results_folder):
        args_str += " --output {} --overwrite".format(results_folder)
        args = args_str.split()
        results_arr = csm_run(args)
        return results_arr


class TestBasic(RunThings):
    test_dir = r"C:\Users\devora\Sources\csm\csm\python\tests\argument_tests\files_for_tests"
    os.chdir(test_dir)
    results_folder = "csm_tests"
    #try:
    #    shutil.rmtree(results_folder)
    #except FileNotFoundError:
    #    pass
    #os.mkdir(results_folder)

    def run_args(self, args_str):
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
        # --remove-hy removes hy. can't be used with select-atoms

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
            file.readline()  # skip the openbabel line, which changes every time
            output_str = file.read()
        with open(os.path.join("expected", "removehyexpected.mol"), 'r') as file:
            file.readline()
            file.readline()  # skip the openbabel line, which changes every time
            expected_str = file.read()

        assert expected_str == output_str

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
            file.readline()  # skip the openbabel line, which changes every time
            output_str = file.read()
        with open(os.path.join("expected", "selectatomsexpected.mol"), 'r') as file:
            file.readline()
            file.readline()  # skip the openbabel line, which changes every time
            expected_str = file.read()

        assert expected_str == output_str

    def test_select_atoms_remove_hy(self):
        # --select-atoms removes specific atoms.
        # --remove-hy removes 'H' atoms.

        cmd = "exact c2 --input 4-helicene.mol --select-atoms 15-19,1,2 --remove-hy"
        results = self.run_args(cmd)
        assert len(results[0][0].molecule) == 7

        # test output
        with open(os.path.join(self.results_folder, "resulting_symmetric_coordinates.mol"), 'r') as file:
            file.readline()
            file.readline()  # skip the openbabel line, which changes every time
            output_str = file.read()
        with open(os.path.join("expected", "selectatomsexpected.mol"), 'r') as file:
            file.readline()
            file.readline()  # skip the openbabel line, which changes every time
            expected_str = file.read()

        assert expected_str == output_str

    def test_ignore_atoms(self):
        # --ignore-atoms removes specific atoms.
        # The test checks that the result is identical to the result by use --select-atoms

        # cmd = "exact c2 --input 4-helicene.mol --select-atoms 15-19,1,2"  # len(_atoms) = 30
        cmd = "exact c2 --input 4-helicene.mol --ignore-atoms 3-14,20-30"
        results = self.run_args(cmd)
        assert len(results[0][0].molecule) == 7

        # test output
        with open(os.path.join(self.results_folder, "resulting_symmetric_coordinates.mol"), 'r') as file:
            file.readline()
            file.readline()  # skip the openbabel line, which changes every time
            output_str = file.read()
        with open(os.path.join("expected", "selectatomsexpected.mol"), 'r') as file:
            file.readline()
            file.readline()  # skip the openbabel line, which changes every time
            expected_str = file.read()

        assert expected_str == output_str

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

        # TODO: output tests

    # output
    def test_legacy(self):
        #why is it printing filename instead of index all of a sudden?
        # --legacy-output prints old style csm format
        cmd = "exact c2 --input squarate.xyz --select-mols 1 --legacy-output --output {}/legacy.txt".format(
            self.results_folder)
        csm_run(cmd.split())
        with open(os.path.join(self.results_folder, "legacy.txt"), 'r') as ofile:
            out = ofile.read()
        with open(os.path.join("expected", "legacy.txt"), 'r') as ofile:
            exp = ofile.read()
        assert out == exp

    def test_json_output(self):
        #same problem as legacy
        # --json-output. only works with --legacy-output
        cmd = "exact c2 --input squarate.xyz --select-mols 1 --legacy-output --output {}/legacy.json --json-output".format(
            self.results_folder)
        csm_run(cmd.split())
        with open(os.path.join(self.results_folder, "legacy.json"), 'r') as ofile:
            out = ofile.read()
        with open(os.path.join("expected", "legacy.json"), 'r') as ofile:
            exp = ofile.read()
        assert out == exp

    def test_print_denorm(self):
        # --print-denorm prints denormalized coordinates for original molecule instead of normalized
        cmd = "exact c2 --input squarate.xyz --select-mols 1 --print-denorm"
        self.run_args(cmd)
        with open(os.path.join(self.results_folder, "initial_normalized_coordinates.xyz"), "r") as file:
            output2 = file.read()

        cmd = "exact c2 --input squarate.xyz --select-mols 1"
        self.run_args(cmd)
        with open(os.path.join(self.results_folder, "initial_denormalized_coordinates.xyz"), "r") as file:
            output1 = file.read()

        assert output1 != output2

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
        assert out == exp

    def test_verbose(self):
        # verbose-- approx only
        cmd = "approx c2 --input squarate.xyz --verbose"
        self.run_args(cmd)
        output_path = os.path.join(self.results_folder, "approx")
        assert os.path.isdir(output_path)
        # todo: because the output contains a variable runtime, running a comparison is a bit tedious, leaving it for now

    # shared stuff

    def test_sn_max(self):
        # --sn-max (relevant only for chirality)
        cmd = "exact ch --input bis(dth)copper(I).mol"
        result1 = self.run_args(cmd)
        cmd = "exact ch --input bis(dth)copper(I).mol --sn-max 2"
        result2 = self.run_args(cmd)
        assert result1[0][0].op_type == 'SN' and result1[0][0].op_order == 4
        assert result2[0][0].op_type == 'CS' and result2[0][0].op_order == 2

    def test_normalize(self):
        cmd = "approx c3 --fibonacci 1 --read-fragments --input reading-fragments\water-6.mol --normalize 0 1 2 3 4 5 6"
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

    def test_keep_structure(self):
        cmd = "exact cs --input 4-helicene.mol --keep-structure"
        results = self.run_args(cmd)
        assert results[0][0].csm == pytest.approx(0.793551, abs=1e-5)

    def test_output_perms(self):
        cmd = "exact cs --input 4-helicene.mol --keep-structure --output-perms"
        self.run_args(cmd)
        out_rows = []
        with open(os.path.join(self.results_folder, "perms.csv"), 'r') as file:
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
    def test_use_chains(self):
        cmd = "approx c3 --input 2rla-s3.pdb --use-chains"
        results = self.run_args(cmd)
        assert len(results[0][0].molecule.chains) == 3
        cmd = "approx c3 --input 2rla-s3.pdb"
        results = self.run_args(cmd)
        assert len(results[0][0].molecule.chains) == 1

    def test_selective(self):
        cmd = "approx c3 --input 2rla-s3.pdb --fibonacci 50 --selective 3"
        results = self.run_args(cmd)

    def test_permute_chains(self):
        cmd = "trivial c3 --input 2rla-s3.pdb"
        results1 = self.run_args(cmd)
        assert results1[0][0].csm == pytest.approx(49.078986, abs=1e-5)

        cmd = "trivial c3 --input 2rla-s3.pdb --permute-chains"
        results2 = self.run_args(cmd)
        assert results2[0][0].csm == pytest.approx(0.009663, abs=1e-5)

    def test_directions(self):
        # TODO: refine detect-outliers test
        # default
        cmd = "approx c3 --input 3alb-gkt4-h.pdb"
        results = self.run_args(cmd)
        assert len(results[0][0].ongoing_statistics['approx']) == 9
        # detect-outliers
        cmd = "approx c3 --input 3alb-gkt4-h.pdb --detect-outliers"
        results = self.run_args(cmd)
        assert len(results[0][0].ongoing_statistics['approx']) == 18
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

    def test_algorithms(self):
        # default
        cmd = "approx c3 --input 2rla-s3.pdb"
        results = self.run_args(cmd)
        assert round(results[0][0].csm, 6) == 0.009663
        stats1 = tuple(results[0][0].ongoing_statistics['approx'][3]['stats']['dirs'][0])
        # greedy
        cmd = "approx c3 --input 2rla-s3.pdb --greedy"
        results = self.run_args(cmd)
        assert round(results[0][0].csm, 6) == 0.009663
        stats2 = tuple(results[0][0].ongoing_statistics['approx'][3]['stats']['dirs'][0])
        # many-chains
        cmd = "approx c3 --input 2rla-s3.pdb --many-chains"
        results = self.run_args(cmd)
        assert round(results[0][0].csm, 6) == 0.009663
        stats3 = tuple(results[0][0].ongoing_statistics['approx'][3]['stats']['dirs'][0])
        # keep-structure
        cmd = "approx c3 --input 2rla-s3.pdb --keep-structure"
        results = self.run_args(cmd)
        assert round(results[0][0].csm, 6) == 0.009663
        stats4 = tuple(results[0][0].ongoing_statistics['approx'][3]['stats']['dirs'][0])

        test = set([stats1, stats2, stats3, stats4])
        assert len(test) == 4

    #parallel
    def test_parallel_mols_in_file(self):
        #TODO
        assert False

class TestFragments(RunThings):
    test_dir = r"C:\Users\devora\Sources\csm\csm\python\tests\argument_tests\files_for_tests"
    os.chdir(test_dir)
    results_folder = "csm_tests"
    #shutil.rmtree(results_folder)
    #os.mkdir(results_folder)

    def run_args(self, args_str):
        return super()._run_args(args_str, self.results_folder)

    def test_read_fragments_baseline(self):
        # --read-fragments reads from .mol, .pdb as chains
        filenames = ["model-endmdl-withIDS.pdb",
                     "model-endmdl-withoutIDS.pdb",
                     "water-6.mol",
                     "water-6.pdb"]

        for filename in filenames:
            cmd = "approx c3 --fibonacci 1 --input reading-fragments"
            cmd += "\\" + filename
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
            cmd = "approx c3 --fibonacci 1 --read-fragments --input reading-fragments"
            cmd += "\\" + filename
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
            cmd = "approx c3 --fibonacci 1 --read-fragments --use-chains --input reading-fragments"
            cmd += "\\" + filename
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
        cmd = "trivial ch --input 2rla-s3.pdb"
        result1 = self.run_args(cmd)
        assert result1[0][0].csm == pytest.approx(0.111737, abs=1e-5)
        assert result1[0][0].overall_statistics["best chirality"] == "CS"

    def test_trivial_4(self):
        cmd = "trivial ch --input 2rla-s3.pdb --permute-chains"
        result2 = self.run_args(cmd)
        assert result2[0][0].csm == pytest.approx(0.111737, abs=1e-5)
        assert result2[0][0].overall_statistics["best chirality"] == "CS"

    def test_exact(self):
        # --sn-max (relevant only for chirality)
        cmd = "exact ch --input bis(dth)copper(I).mol" #matches
        result1 = self.run_args(cmd)
        assert result1[0][0].csm == pytest.approx(0, abs=1e-5)
        assert result1[0][0].overall_statistics["best chirality"] == "S4"

    def test_exact_2(self):
        cmd = "exact ch --input bis(dth)copper(I).mol --sn-max 2" #matches
        result2 = self.run_args(cmd)
        assert result2[0][0].csm == pytest.approx(4.394381, abs=1e-5)
        assert result2[0][0].overall_statistics["best chirality"] == "CS"

    def test_approx(self):
        cmd = "approx ch --input ferrocene.xyz --babel-bond" #matches
        results = self.run_args(cmd)
        assert results[0][0].csm == pytest.approx(0.241361, abs=1e-5)
        assert results[0][0].overall_statistics["best chirality"] == "CS"

    def test_approx_2(self):
        cmd = "approx ch --input bis(dth)copper(I).mol" #matches
        result1 = self.run_args(cmd)
        assert result1[0][0].csm == pytest.approx(0, abs=1e-5)
        assert result1[0][0].overall_statistics["best chirality"] == "S4"

    def test_approx_3(self):
        cmd = "approx ch --input bis(dth)copper(I).mol --sn-max 2" #matches
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
