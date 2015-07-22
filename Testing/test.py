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

    """def test_12(self):
        self.one_test('test12')"""

    """ Test 13 gives wrong result in C++ as a result of different versions of OpenBabel. We can't use it.
    def test_13(self):
        self.one_test('test13')
    """

    def test_16(self):
        self.one_test('test16')

    def test_17(self):
        self.one_test('test17')


class LongCSMCheckOutput(TestCase):
    def one_test(self, test_dir):
        self.assertTrue(run_test_python(os.path.join(config.TEST_PATH, test_dir)))


    def test_14(self):
        self.one_test('test14')

    """Test 15 is takes a lot of time"""
    #def test_15(self):
    #    self.one_test('test15')

    """def test_18(self):
        self.one_test('test18')"""

    """Test 19 is takes a lot of time"""
    #def test_19(self):
    #    self.one_test('test19')
