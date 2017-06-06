def pytest_addoption(parser):
    parser.addoption("--dir", action="append", default=[],
        help="list of test directories to pass to test functions")
    parser.addoption("--parent", action="append", default=[],
                     help="list of parent directories, whose child directories will be passed to test functions")

def pytest_generate_tests(metafunc):
    pass