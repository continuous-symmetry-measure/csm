import os
import json
from tests.utils.run_test import Expected

def pytest_addoption(parser):
    parser.addoption("--dir", action="append", default=[],
        help="list of test directories to pass to test functions")
    parser.addoption("--parent", action="append", default=[],
                     help="list of parent directories, whose child directories will be passed to test functions")

def pytest_generate_tests(metafunc):
    dirs=metafunc.config.option.dir[:]
    for test_folder in metafunc.config.option.parent:
        for folder in os.listdir(test_folder):
            mydir=os.path.join(test_folder,folder)
            if os.path.isdir(mydir):
                dirs.append(mydir)
    params=get_run_tuples(dirs)


    if 'args' in metafunc.fixturenames:
        metafunc.parametrize("args, molecule, expected, equiv", params)



def get_run_tuples(dirs):
    params=[]
    for mydir in dirs:
        inputjson=os.path.join(mydir, 'input.json')
        outputjson=os.path.join(mydir, 'output.json')
        molecule_folder = os.path.join(mydir, 'molecules')



        with open(inputjson) as f:
            in_dict=json.load(f)
        with open(outputjson) as f:
            output_dict = json.load(f)


        for key in in_dict["runs"]:
            for molecule in os.listdir(molecule_folder):
                molfile=os.path.join(molecule_folder, molecule)
                mol_index = molecule.split(".")[0]
                e = output_dict[key][mol_index]
                expected = Expected(e)
                my_tuple=(in_dict["runs"][key], molfile, expected, in_dict['equiv_perms'])
                params.append(my_tuple)

    return params

