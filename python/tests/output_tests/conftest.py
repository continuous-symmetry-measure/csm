import os
import shutil

version="0.20.8"
session_folder=r"C:\Users\devora\Sources\temp\csm_tests"
test_folder=r"C:\Users\devora\Sources\csm\csm\python\tests\output_tests\files_for_tests"

def pytest_generate_tests(metafunc):
    #clear existing result files, because otherwise you end up with timestamped results not being checked
    #edit: not necessary because i make timestamps...
    #try:
    #    shutil.rmtree(session_folder)
    #except FileNotFoundError:
    #    pass
    #os.makedirs(session_folder)

    #parametrize
    if 'test_name' in metafunc.fixturenames:
        metafunc.parametrize("test_name", [dir for dir in os.listdir(test_folder)])