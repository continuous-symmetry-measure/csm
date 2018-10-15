import os

version="0.20.8"
session_folder=r"C:\Users\devora\Sources\temp\csm_tests"
test_folder=r"C:\Users\devora\Sources\csm\tests"

def pytest_generate_tests(metafunc):
    if 'test_name' in metafunc.fixturenames:
        metafunc.parametrize("test_name", [dir for dir in os.listdir(test_folder)])