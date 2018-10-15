import os

import datetime
import pytest

from csm.main.csm_run import csm_run
from output_tests.conftest import session_folder, version, test_folder

def standard_folder(test_name):
    return os.path.join(test_folder, test_name, "test_standard")

def close_enough(x, y):
    return x==y

def check_correct(result, test_standard):
    try:
        assert close_enough(float(result), float(test_standard))
    except ValueError:
        assert result == test_standard

def check_file(result_file, test_standard_file):
    with open(result_file, 'r') as r, open(test_standard_file, 'r') as s:
        for line_r, line_s in zip(r,s):
            split_r=line_r.strip().split()
            split_s=line_s.strip().split()
            for res, stand in zip (split_r, split_s):
                check_correct(res, stand)


class CheckFolder:
    timestamp=str(datetime.datetime.now().timestamp())[-11:].replace(".", "")

    @pytest.fixture(scope="class")
    def test_name(self):
        return

    def result_folder(self, test_name):
        return os.path.join(session_folder, test_name, self.timestamp)

    def test_run(self, test_name):
        pass

    def call_checker(self, filename, test_name):
        check_file(os.path.join(self.result_folder(test_name), filename),
                   os.path.join(standard_folder(test_name), filename))

    def test_files(self, test_name):
        standard=standard_folder(test_name)
        result=self.result_folder(test_name)
        for root, dirs, files in os.walk(standard):
            for file in files:
                if file=="extra.txt":
                    continue
                self.call_checker(file, test_name)
                print(file, "passed")
            break
        for dir in dirs:
            for file in os.listdir(os.path.join(standard, dir)):
                check_file(os.path.join(result, dir, file),
                           os.path.join(standard, dir, file))
                print(file, "passed")




class xTest:
    def test_csm(self, test_name):
        filename="csm.txt"
        self.call_checker(filename, test_name)

    def test_directional(self, test_name):
        filename="directional.txt"
        self.call_checker(filename, test_name)

    def test_extra(self, test_name):
        filename="extra.txt"
        self.call_checker(filename, test_name)

    def test_permutation(self, test_name):
        filename="permutation.txt"
        self.call_checker(filename, test_name)

    def test_version(self, test_name):
        filename="version.txt"
        assert "CSM VERSION: "+version

    def test_old_csm_output(self, test_name):
        foldername="old-csm-output"
        res_fold=os.path.join(self.result_folder(test_name), foldername)
        stand_fold=os.path.join(standard_folder(test_name), foldername)
        files=os.listdir(stand_fold)
        for file in files:
            check_file(os.path.join(res_fold, file),
                       os.path.join(stand_fold, file))


    def test_initial_coordinates(self, test_name):
        pass

    def test_symmetric_coordinates(self, test_name):
        pass
