from unittest import TestCase

__author__ = 'zmbq'

from helpers import run_test
import config
import os.path

class CheckCSMOutput(TestCase):
    def one_test(self, test_dir):
        self.assertTrue(run_test(os.path.join(config.TEST_PATH, test_dir)))

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



if __name__=='__main__':
    run_test(os.path.join(config.TEST_PATH, 'test2'))

