import os
import json
from tests.utils.run_test import Expected

def pytest_addoption(parser):
    parser.addoption("--dir", type=str)

def pytest_generate_tests(metafunc):
    mydir=metafunc.config.option.dir
    inputjson=os.path.join(mydir, 'input.json')
    outputjson=os.path.join(mydir, 'output.json')
    molecule_folder = os.path.join(mydir, 'molecules')

    params=[]

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
            my_tuple=(in_dict["runs"][key], molfile, expected)
            params.append(my_tuple)


    metafunc.parametrize("args, molecule, expected", params)
