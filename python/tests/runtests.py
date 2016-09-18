import sys
import os
from tests.utils.run_test import run_test


if __name__ == '__main__':
    args = sys.argv[1:]
    test_folder=args[0]
    for folder in os.listdir(test_folder):
        try:
            run_test(os.path.join(test_folder, folder))
        except NotADirectoryError:
            pass  # we dont care