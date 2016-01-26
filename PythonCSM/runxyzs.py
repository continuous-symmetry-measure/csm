from molecule.molecule import Molecule
from openbabel import OBAtomAtomIter, OBConversion, OBMol
from molecule.atom import GetAtomicSymbol, Atom
from a_calculations.csm_calculations import exact_calculation
import numpy as np
import csv

def split(filename):
    mol_dict={}
    with open(filename, 'r') as file:
        chunks=file.read().split("MOL_INDEX=")
        for item in chunks:
            try:
                    firstcoords=item.find("INITIAL STRUCTURE COORDINATES")
                    rescoords=item.find("RESULTING STRUCTURE COORDINATES")
                    resparams=item.find("DIRECTIONAL COSINES")

                    settings=item[:firstcoords]
                    input=item[firstcoords+30:rescoords]
                    output=item[rescoords+32:resparams]
                    results=item[resparams:]
                    (index, symmetryline, scalingfactor)=settings.strip().split("\n")
                    symm=symmetryline[:2]
                    csm=symmetryline[symmetryline.find(":")+1:]
                    (filler, blank, dir, blank2, filler2, blank3, perm) = results.strip().split("\n")
                    mol_dict[index,symm]={"in_coord":input, "out_atoms":output, "csm":csm, "scalingfactor":scalingfactor, "dir":dir, "perm":perm}
            except:
                continue
    return mol_dict




def xyz_split(filename):
    '''
    :param filename: name of xyz file containing multiple molecules, to be split
    :return: dictionary owith key: index of molecule and val: molecule string (params are not included)
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

    #fetching the result and validation data:
    r_perm=[p + 1 for p in result.perm]
    v_perm= [int(i) for i in validate["perm"].split()]
    r_csm=result.csm
    v_csm=float(validate["csm"])
    r_dir= result.dir
    v_dir= [float(i) for i in validate["dir"].split()]
    r_atoms=result.outAtoms
    v_mol=Molecule.from_string(validate["out_atoms"], "xyz", False)
    v_atoms= np.asarray([np.asarray(atom.pos) for atom in v_mol.atoms])

    #build the dictionary:
    res={}
    res["molecule id"]=mol_index
    res["symmetry"]=symm
    res["csm result"]=r_csm
    res["dir result"]=r_dir
    res["perm result"]=r_perm
    res["atoms result"]=r_atoms
    res["csm expected"]=v_csm
    res["dir expected"]=v_dir
    res["perm expected"]=v_perm
    res["atoms expected"]=v_atoms


    #check the status:
    failed=False
    res["status"]="Failed"
    res["message"]=""
    #check csm:
    if abs(r_csm -v_csm)>.0001: #greater than 4 decimal points: the validation stops, hence the result will always be bigger (since it continue for several more decimal points
        failed=True
        res["message"]+= "CSM does not match "

    #chck symmetric structure
    atom_pos_not_match=[]
    for i in range(len(r_atoms)):
        for index in range(3):
            val1= r_atoms[i][index]
            val2= v_atoms[i][index]
            if abs(val1 -val2)>.01:
                atom_pos_not_match.append((i,index,val1,val2))
    if atom_pos_not_match:
        failed=True
        res["message"]+= "Symmetric structure does not match"

    if failed:
        return res

    res["status"]="OK"
    if r_perm != v_perm:
        res["message"]+="different perm "

    dir_not_match=False
    dir_not_times_minus=False
    for i in range(3):
        val1= r_dir[i]
        val2= (v_dir[i])
        if abs(val1 -val2)>.0001:
            dir_not_match=True
            if abs(-val1 -val2)>.0001:
                dir_not_times_minus=True
    if dir_not_match:
        if dir_not_times_minus:
            res["message"]+="unexpected dir"
        else:
            res["message"]+="dir reversed sign"

    return res




def runtest(molecule_file, symmetry_file, result_file):
    molecules = xyz_split(molecule_file)
    validation=split(result_file)
    with open(r'C:\Users\dev\Documents\Chelem\csm\PythonCSM\testing.csv', 'w') as csvfile:
        fieldnames = ['molecule id','symmetry','status','message','csm result','csm expected','atoms result','atoms expected','dir result','perm result','dir expected','perm expected']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for index, symm in validation:
            try:
                symmetry=symm.lower()
                xyz = molecules[int(index)]
                molecule = Molecule.from_string(xyz, "xyz")
                result = exact_calculation(symmetry, molecule)
                test=result.operationType
                validate=validation[index, symm]
                res=check_result(result, validate, index, symm)
                writer.writerow(res)
            except:
                pass


def run():
    directory=r'C:\Users\dev\Documents\Chelem\csm\mols_for_Itay'
    molfile = directory+ r'\input\input_1_methane_csm\methane-test1.xyz'
    symmfile = directory+ r'\expected_output\expected_output_1_methane_csm\sym.txt'
    resfile = directory + r'\expected_output\expected_output_1_methane_csm\csmresults.log'
    runtest(molfile, symmfile, resfile)


run()
