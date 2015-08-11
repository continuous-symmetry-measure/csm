from unittest import TestCase

__author__ = 'zmbq'

from helpers import run_test_python
import config
import os.path

class CheckCSMOutput(TestCase):
    def one_test(self, test_dir):
        self.assertTrue(run_test_python(os.path.join(config.TEST_PATH, test_dir)))

    def test_1(self):
        self.one_test('test1')

    def test_2(self):
        self.one_test('test2')

    def test_3(self):
        self.one_test('test3')

    def test_4(self):
        self.one_test('test4')

    # Test 5 was never created...

    def test_6(self):
        self.one_test('test6')

    def test_7(self):
        self.one_test('test7')

    def test_direction(self):
        self.one_test('test_direction')

    def test_permutation(self):
        self.one_test('test_permutation')

    def test_8(self):
        self.one_test('test8')

    def test_9(self):
        self.one_test('test9')

    """Open babel can not open the molecule file for test 12
    def test_12(self):
        self.one_test('test12')"""

    def test_15(self):
        self.one_test('test15')

    def test_16(self):
        self.one_test('test16')

    def test_17(self):
        self.one_test('test17')

    def test_20(self):
        self.one_test('test20')

    def test_21(self):
        self.one_test('test21')

    # in test 22 the resulting permutation and direction might be different (?) - see output files in test22 folder
    def test_22(self):
        self.one_test('test22')


class LongCSMCheckOutput(TestCase):
    def one_test(self, test_dir):
        self.assertTrue(run_test_python(os.path.join(config.TEST_PATH, test_dir)))

    def test_13(self):
        self.one_test('test13csm')

    def test_14(self):
        self.one_test('test14')

    """def test_18(self):
        self.one_test('test18')"""

    """Test 19 is takes a lot of time"""
    #def test_19(self):
    #    self.one_test('test19')
