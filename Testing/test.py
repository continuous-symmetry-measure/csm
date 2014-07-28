from unittest import TestCase

__author__ = 'zmbq'

from helpers import run_test
import config
import os.path

class CheckTestInstances(TestCase):
    def test_entire_folder(self):
        for test_dir in os.listdir(config.TEST_PATH):
            self.assertTrue(run_test(os.path.join(config.TEST_PATH, test_dir)))

if __name__=='__main__':
    run_test(os.path.join(config.TEST_PATH, 'test2'))

