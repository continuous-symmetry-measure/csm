from calculations.permuters import is_legal_perm
from collections import OrderedDict
from input_output.arguments import get_operation_data

from molecule.molecule import Molecule
from openbabel import OBAtomAtomIter, OBConversion, OBMol
from molecule.atom import GetAtomicSymbol, Atom
from calculations.csm_calculations import exact_calculation
import numpy as np
import csv
import re
import logging
import os.path
import glob


opcode_data = {
    "cs": ('CS', 2, "MIRROR SYMMETRY"),
    "ci": ('CI', 2, "INVERSION (S2)"),
    "ch": ('CH', 2, "CHIRALITY"),
    "c2": ('CN', 2, "C2 SYMMETRY"),
    'c3': ('CN', 3, "C3 SYMMETRY"),
    'c4': ('CN', 4, "C4 SYMMETRY"),
    'c5': ('CN', 5, "C5 SYMMETRY"),
    'c6': ('CN', 6, "C6 SYMMETRY"),
    'c7': ('CN', 7, "C7 SYMMETRY"),
    'c8': ('CN', 8, "C8 SYMMETRY"),
    's2': ('SN', 2, "S2 SYMMETRY"),
    's4': ('SN', 4, "S4 SYMMETRY"),
    's6': ('SN', 6, "S6 SYMMETRY"),
    's8': ('SN', 8, "S8 SYMMETRY")
}

def split(filename):
    mol_dict = OrderedDict()
    with open(filename, 'r') as file:
        chunks = file.read().split("MOL_INDEX=")
        for item in chunks:
            try:
                firstcoords = item.find("INITIAL STRUCTURE COORDINATES")
                rescoords = item.find("RESULTING STRUCTURE COORDINATES")
                resparams = item.find("DIRECTIONAL COSINES")

                settings = item[:firstcoords]
                input = item[firstcoords + 30:rescoords]
                output = item[rescoords + 32:resparams]
                results = item[resparams:]
                (index, symmetryline, scalingfactor) = settings.strip().split("\n")
                # indexcheck=re.search('[a-zA-Z]+', indexline)
                # index=indexline[:indexcheck.start()].strip()
                symm = symmetryline[:2]
                csm = symmetryline[symmetryline.find(":") + 1:]
                (filler, blank, dir, blank2, filler2, blank3, perm) = results.strip().split("\n")
                mol_dict[index, symm] = {"in_coord": input, "out_atoms": output, "csm": csm,
                                         "scalingfactor": scalingfactor, "dir": dir, "perm": perm}
            except:
                continue
    return mol_dict


def xyz_split(filename):
    '''
    :param filename: name of xyz file containing multiple molecules, to be split
    :return: dictionary with key: index of molecule and val: molecule string (params are not included)
    '''
    mol_dict = {}
    index = 0
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
                mol_dict[index] = mol
            except:
                continue
    return mol_dict


def get_symmetries_list(filename):
    with open(filename, 'r') as file:
        args = []
        for line in file:
            CSM = line.find("csm") + 3
            IN = line.find("__IN__")
            arg = line[CSM:IN]
            cleaned = arg.replace(" ", "")
            args.append(cleaned)
    return args


def check_result(result, validate, mol_index, symm):
    # fetching the result and validation data:
    r_perm = [p + 1 for p in result.perm]
    v_perm = [int(i) for i in validate["perm"].split()]
    r_csm = result.csm
    v_csm = float(validate["csm"])
    r_dir = result.dir
    v_dir = [float(i) for i in validate["dir"].split()]
    r_atoms = result.symmetric_structure
    v_mol = Molecule.from_string(validate["out_atoms"], "xyz", False)
    v_atoms = np.asarray([np.asarray(atom.pos) for atom in v_mol.atoms])

    # build the dictionary:
    res = {}
    res["molecule id"] = mol_index
    res["symmetry"] = symm
    res["csm result"] = r_csm
    res["dir result"] = r_dir
    res["perm result"] = r_perm
    res["atoms result"] = r_atoms
    res["csm expected"] = v_csm
    res["dir expected"] = v_dir
    res["perm expected"] = v_perm
    res["atoms expected"] = v_atoms

    # check the status:
    failed = False
    res["status"] = "Failed"
    res["message"] = []

    # check csm:
    if abs(r_csm - v_csm) > .0001:  # greater than 4 decimal points: the validation stops, hence the result will always be bigger (since it continue for several more decimal points
        failed = True
        res["message"].append("CSM mismatch")

    # chck symmetric structure
    atom_pos_not_match = []
    for i in range(len(r_atoms)):
        for index in range(3):
            val1 = r_atoms[i][index]
            val2 = v_atoms[i][index]
            if abs(val1 - val2) > .01:
                atom_pos_not_match.append((i, index, val1, val2))
    if atom_pos_not_match:
        failed = True
        res["message"].append("Structure mismatch")

    res['verify']=[]
    # check perm
    if r_perm != v_perm:
        res['verify']=True
        res["message"].append("Perm mismatch")

    # check dir
    dir_not_match = False
    dir_not_times_minus = False
    for i in range(3):
        val1 = r_dir[i]
        val2 = (v_dir[i])
        if abs(val1 - val2) > .0001:
            dir_not_match = True
            if abs(-val1 - val2) > .0001:
                dir_not_times_minus = True
    if dir_not_match:
        if dir_not_times_minus:
            res["message"].append("Dir mismatch")
        else:
            res["message"].append("Dir negative")

    if failed:
        return atom_pos_not_match, res
    res["status"] = "OK"
    return atom_pos_not_match, res


def runtests(molecule_file, symmetry_file, result_file, directory, name):
    print(name)
    molecules = xyz_split(molecule_file)
    validation = split(result_file)
    filename = os.path.join(directory, name + ".csv")
    with open(filename, 'w') as csvfile:
        fieldnames = ['molecule id', 'symmetry', 'status', 'message', 'verify', 'csm result', 'csm expected', 'atoms result',
                      'atoms expected', 'perm result', 'perm expected', 'dir result', 'dir expected']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, lineterminator='\n')
        writer.writeheader()
        for index, symm in validation:
            print (index, symm, '...', end='')
            oursymm = 'CS' if symm=='MI' else symm  # CS is marked as MI (Mirror)
            operation = get_operation_data(oursymm)
            xyz = molecules[int(index)]
            molecule = Molecule.from_string(xyz, "xyz")
            validate = validation[index, symm]
            result = exact_calculation(operation.type, operation.order, molecule)
            atom_mismatch, res = check_result(result, validate, index, symm)
            if res['verify']:
                molecule2 = Molecule.from_string(xyz, "xyz")
                expected_perm=[p-1 for p in res['perm expected']]
                verify = exact_calculation(operation.type, operation.order, molecule2,expected_perm)
                atom_mismatch2, res2 = check_result(verify, validate, index, symm)
                if res2['status']!= 'OK':
                    res['verify']=(res2['message'], res2['csm result'], res2['atoms result'])
                else:
                    res['verify']="Verified"
                if not is_legal_perm(expected_perm, molecule2):
                    res['verify']+="illegal original perm"
            writer.writerow(res)
            print()
    print("done")


def test_individuals(molfile, symmfile, permfile, resfile):

    for filename in glob.glob('*.txt'):
        result = exact_calculation(symmetry, molecule, perm=perm)



def run():
    # test_individual()

    # Init logging
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%a, %d %b %Y %H:%M:%S',
                        filename='runxyz.log',
                        filemode='w')
    #directory=r'C:\Users\dev\Documents\Chelem\csm'
    directory = r'../test_cases/inbal'

    name = "methane_test"
    molfile = os.path.join(directory, 'input/input_1_methane_csm/methane-test1.xyz')
    symmfile = os.path.join(directory, 'expected_output/expected_output_1_methane_csm/sym.txt')
    resfile = os.path.join(directory,  'expected_output/expected_output_1_methane_csm/csmresults.log')
    #test_individuals(molfile, symmfile, resfile, directory, name)

    name = "biphenyl_test"
    molfile = os.path.join(directory, r'/input/input_2_biphenyl/biphenyl_test.xyz')
    symmfile = os.path.join(directory, r'/expected_output/expected_output_2_biphenyls/sym.txt')
    resfile = os.path.join(directory, r'/expected_output/expected_output_2_biphenyls/csmresults.log')
    #runtests(molfile, symmfile, resfile, directory, name)

    name = "cyclopentadiene_test"
    molfile = os.path.join(directory, r'input/input_3_cyclopentadiene/cyclopentadiene-test.xyz')
    symmfile = os.path.join(directory, r'expected_output/expected_output_3_cyclopentadiene/sym.txt')
    resfile = os.path.join(directory, r'expected_output/expected_output_3_cyclopentadiene/csmresults.log')
    runtests(molfile, symmfile, resfile, directory, name)

    name = "4cluster_test"
    molfile = os.path.join(directory, r'/input/input_4-clusters/W_Au12_optimized_B3P86.xyz')
    symmfile = os.path.join(directory, r'/expected_output/expected_output_4_clusters/sym.txt')
    resfile = os.path.join(directory, r'/expected_output/expected_output_4_clusters/csmresults.log')
    #runtests(molfile, symmfile, resfile, directory, name)


if __name__ == '__main__':
    run()
