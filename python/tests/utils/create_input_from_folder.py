import sys
import os
import glob
import json

test_cases_path=r'D:\UserData\devora\Sources\csm\python\tests\test_cases'
sagiv_path=r'D:\UserData\devora\Sources\csm\test_cases\sagiv\test_cases\clean_test_cases'

def get_file_names(dir):
    if not os.path.exists(dir):
        raise NotADirectoryError
    argfile=os.path.join(dir, 'sym.txt')
    permfile=os.path.join(dir, 'eq-perm.txt')
    molfiles=[]
    for molfilename in glob.glob(os.path.join(dir, '*.xyz')):
        molfiles.append(molfilename)
    if len(molfiles)>1:
        print("more than one zyx file in input folder?")
    molfile=molfiles[0]
    return argfile, permfile, molfile

def xyz_split(filename, foldername):
    '''
    :param filename: name of file containing xyz molecules to be split
    :param foldername: name of folder to save split xyz molecules in
    :return:
    '''
    index = 0
    if not os.path.exists(foldername):
        os.makedirs(foldername)

    with open(filename, 'r') as file:
        for line in file:
            try:
                mol = ""
                num = int(line)
                mol = line
                mol += file.readline()
                for x in range(num):
                    mol += file.readline()
                index += 1
                filepath=os.path.join(foldername, str(index)+".xyz")
                with open(filepath, 'w') as f:
                    f.write(mol)

            except:
                continue


def get_equiv_perms(permfile):
    perms_dict={}
    try:
        with open(permfile) as f:
            test = f.readline()
            symm=test.strip().split("\t")[2]
            perms_dict[symm]=[]
            while True:
                key=f.readline().strip() #[int(x) for x in f.readline().strip().split()]
                if not key:
                    break
                key_arr=[int(x)-1 for x in key.strip().split()]
                perms_dict[symm].append([key_arr])
                while True:
                    test=f.readline().strip()
                    if not test: #reached eof
                        break
                    if test[0]=='p':
                        symm = test.split("\t")[2]
                        if symm not in perms_dict:
                            perms_dict[symm] = []
                        break
                    key_arr = [int(x) - 1 for x in test.split()]
                    perms_dict[symm][-1].append(key_arr)#[int(x) for x in test.strip().split()])

    except FileNotFoundError:
        pass #we don't care, eq perms not required

    finally:
        return perms_dict

def create_json(dir, test_name):
    argfile, permfile, molfile= get_file_names(dir)

    print("-----------")
    runs_dict={}

    equiv=get_equiv_perms(permfile)

    with open(argfile) as line_args:
        for index, line in enumerate(line_args):
            args=line.split()
            index_str=str(index+1)
            if index<10:
                index_str="0"+index_str
            run_name="L"+index_str+"_"+args[1]
            runs_dict[run_name]=args[1:]

    in_dict={'moleculefile':molfile,
             'equiv_perms':equiv,
             'runs':runs_dict
    }

    save_json(test_name, in_dict)


def save_json(test_name, in_dict):
    new_test_folder=os.path.join(test_cases_path,test_name)
    if not os.path.exists(new_test_folder):
        os.makedirs(new_test_folder)
    else:
        print("as a heads up, overwriting existing folder")

    xyz_split(in_dict['moleculefile'], os.path.join(new_test_folder, 'molecules'))

    with open(os.path.join(new_test_folder, 'input.json'), 'w') as jsonfile:
        json.dump(in_dict, jsonfile)




def main_script(testpath):
    input=r'input\input_1'
    in_path=os.path.join(testpath, input)
    test_name=os.path.basename(os.path.normpath(testpath))
    create_json(in_path, test_name)



for folder in os.listdir(sagiv_path):
    try:
        main_script(os.path.join(sagiv_path, folder))
    except NotADirectoryError:
        pass #we dont care

