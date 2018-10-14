import json
import os

import sys

from csm.molecule.molecule import MoleculeReader, Molecule
from csm.input_output.readers import read_from_sys_std_in
from csm.input_output.formatters import csm_log as print

def read_molecules(**kwargs):
    input_name=kwargs["in_file_name"]
    if os.path.isdir(input_name):
        kwargs.pop('in_file_name')
        mols=[]
        files=[]
        for directory, subdirectories, files1 in os.walk(input_name):
            files+=files1
            break

        files.sort()
        for file_name in files:
                if file_name=="sym.txt":
                    continue
                mol_file= os.path.join(input_name, file_name)
                try:
                    mol = MoleculeReader.from_file(mol_file, **kwargs)
                    mols.append(mol)
                except Exception as e:
                    if mol_file[-3:] in ["mol", "pdb", "sdf", "xyz", "csm", "cif"]:
                        print("failed to create a molecule from", file_name, str(e))
                    #else:
                    #    print(file_name, "was not read")

    elif not os.path.isfile(input_name):
        raise ValueError("invalid file/folder name for molecule")
    else: #file
        mols=MoleculeReader.multiple_from_file(**kwargs)
    sys.stderr.flush()

    if kwargs['select_mols']:
        try:
            mols=[mols[i] for i in kwargs['select_mols']]
        except IndexError:
            raise IndexError("You have selected more molecules than you have input")
    return mols

def read_mols_from_std_in():
    raw_json = read_from_sys_std_in()
    less_raw_json = json.loads(raw_json)
    molecules = [Molecule.from_dict(json_dict) for json_dict in less_raw_json]
    return molecules

def read(**dictionary_args):
    mols = read_molecules(**dictionary_args)
    sys.stdout.write(json.dumps([mol.to_dict() for mol in mols], indent=4))
    return mols
