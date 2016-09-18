import sys
import os
import glob
import json
from csm.molecule.molecule import Molecule

test_cases_path=r'D:\UserData\devora\Sources\csm\python\tests\test_cases'
sagiv_path=r'D:\UserData\devora\Sources\csm\test_cases\sagiv\test_cases\clean_test_cases'

def get_file_names(dir):
    if not os.path.exists(dir):
        raise NotADirectoryError
    resfile=os.path.join(dir, 'csmresults.log')
    return resfile


def find_between( s, first, last):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""



def create_json(dir, test_name):
    resfile= get_file_names(dir)

    with open(resfile) as f:
        molecules=f.read().split('MOL_')

    while molecules[0]=='':
        molecules=molecules[1:]

    mol_dict={}

    for mol in molecules:
        try:
            lines=mol.split("\n")
            run_name=lines[0].split()[-1]

            run_name=run_name[run_name.index('=')+1:]
            mol_num=int(find_between(mol, "INDEX=", "SYM"))
            symm= lines[1].split()[0]#find_between(mol, "\n", "SYMMETRY").strip()
            csm= float(lines[1].split()[-1])#float(find_between(mol, "SYMMETRY: ", "SCALING"))
            scale=float(find_between(mol,"FACTOR: ", "INITIAL"))

            normalized_molecule=find_between(mol, "INITIAL STRUCTURE COORDINATES\n", "RESULTING STRUCTURE COORDINATES")
            symmetric_molecule=find_between(mol,"RESULTING STRUCTURE COORDINATES\n", "DIRECTIONAL")

            dir=[float(x) for x in find_between(mol, "DIRECTIONAL COSINES:\n\n",  "\n").strip().split()]
            permutation=mol[mol.index("PERMUTATION:")+len("PERMUTATION:"):].strip('\n').split()
            perm=[int(x) for x in permutation]

            if run_name not in mol_dict:
                mol_dict[run_name]={}
            mol_dict[run_name][mol_num]={
                'csm':csm,
                'dir':dir,
                'perm':perm,
                'symmetric':symmetric_molecule,
                'normalized':normalized_molecule,
            }
        except:
            print(resfile)
            raise

    save_json(test_name, mol_dict)



def save_json(test_name, in_dict):
    new_test_folder=os.path.join(test_cases_path,test_name)
    if not os.path.exists(new_test_folder):
        os.makedirs(new_test_folder)
    else:
        print("as a heads up, overwriting existing folder")

    with open(os.path.join(new_test_folder, 'output.json'), 'w') as jsonfile:
        json.dump(in_dict, jsonfile)

    #refuses to save molecule indexes as ints for some reason
    #with open(os.path.join(new_test_folder, 'output.json')) as jsonfile:
    #    test=json.load(jsonfile)
    #    hi=1






def main_script(testpath):
    addon=r'expected_output\expected_output_1'
    out_path=os.path.join(testpath, addon)
    test_name=os.path.basename(os.path.normpath(testpath))
    create_json(out_path, test_name)



for folder in os.listdir(sagiv_path):
    try:
        main_script(os.path.join(sagiv_path, folder))
    except NotADirectoryError:
        pass #we dont care

