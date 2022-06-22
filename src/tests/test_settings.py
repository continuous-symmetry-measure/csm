import os
test_dir = os.path.join(os.path.dirname(__file__), "argument_tests", "files_for_tests")
test_dir = r"C:\Sources\CSM\csm\src\tests\test-result"
try:
    from local_settings import *
except:
    pass