# Try the openbabel bindings, to debug loading openbabel on github

import os

def print_env(message):
    print(message)
    print('****************')
    paths = os.environ['PATH'].split(';')
    for path in paths:
        print(path)
    print()

print("Trying openbabel bindings")
print_env("Environment at start")
from csm.main.openbabel_fix import prepare_openbabel
prepare_openbabel()
from openbabel.openbabel import OBConversion
print("openbabel.openbabel.OBConversion imported")
conv = OBConversion()
print("PDB format available: ", conv.SetInFormat("pdb"))
print("Done")