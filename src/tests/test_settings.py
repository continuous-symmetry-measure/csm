import os
test_dir = os.path.join(os.path.dirname(__file__), "argument_tests", "files_for_tests")

try:
    from local_settings import *
except:
    pass