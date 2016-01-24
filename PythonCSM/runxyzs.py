from molecule.molecule import Molecule
from openbabel import OBAtomAtomIter, OBConversion, OBMol
from molecule.atom import GetAtomicSymbol, Atom
from a_calculations.csm_calculations import exact_calculation
import numpy as np

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


def check_result(result, validate, index, symm):
    passed=True
    r_perm=[p + 1 for p in result.perm]
    v_perm= [int(i) for i in validate["perm"].split()]
    r_csm=result.csm
    v_csm=float(validate["csm"])

    r_dir= result.dir
    v_dir= [float(i) for i in validate["dir"].split()]
    dir_not_match=False
    dir_not_times_minus=False
    for i in range(3):
        val1= r_dir[i]
        val2= (v_dir[i])
        if abs(val1 -val2)>.0001:
            dir_not_match=True
            if abs(-val1 -val2)>.0001:
                dir_not_times_minus=True

    r_atoms=result.outAtoms
    v_mol=Molecule.from_string(validate["out_atoms"], "xyz")
    v_atoms= np.asarray([np.asarray(atom.pos) for atom in v_mol.atoms])
    atom_pos_not_match=False
    for i in range(len(r_atoms)):
        for index in range(3):
            val1= r_atoms[i][index]
            val2= v_atoms[i][index]
            if abs(val1 -val2)>.0001:
                atom_pos_not_match=True






    if abs(r_csm -v_csm)>.0001: #greater than 4 decimal points: the validation stops, hence the result will always be bigger (since it continue for several more decimal points
        passed=False
        print("failed: molecule", index,"(", symm, ")\t perms:", r_perm, v_perm, "\t csm:", r_csm, v_csm)
        return
    if not atom_pos_not_match:
        print("structure matched")
        return
        #print("symmetric molecule is different:\t validate: ", v_atoms, "\tresult: ", r_atoms)

    passtext="passed"
    if r_perm != v_perm:
        passtext+=", but received different perms than expected: "+str(v_perm)+","+str(r_perm)
    if dir_not_match:
        if dir_not_times_minus:
            passtext+=", but received different dir than expected: "+str(v_dir)+","+str(r_dir)
        else:
            passtext+=", but dir has reversed sign"



    if passed:
        print(passtext)



def run():
    filename = r'C:\Users\devora.witty\Sources\csm\Testing\testing for proteins\tests\mols_for_Itay\input\input_1_methane_csm\methane-test1.xyz'
    molecules = xyz_split(filename)
    filename = r'C:\Users\devora.witty\Sources\csm\Testing\testing for proteins\tests\mols_for_Itay\expected_output\expected_output_1_methane_csm\sym.txt'
    symmetries = get_symmetries_list(filename)
    filename = r'C:\Users\devora.witty\Sources\csm\Testing\testing for proteins\tests\mols_for_Itay\expected_output\expected_output_1_methane_csm\csmresults.log'
    validation=split(filename)
    for index, symm in validation:
        try:
            symmetry=symm.lower()
            xyz = molecules[int(index)]
            molecule = Molecule.from_string(xyz, "xyz")
            result = exact_calculation(symmetry, molecule)
            test=result.operationType
            validate=validation[index, symm]
            check_result(result, validate, index, symm)
        except:
            pass




            # filename=r'C:\Users\devora.witty\Sources\csm\Testing\testing for proteins\tests\mols_for_Itay\input\input_4-clusters\W_Au12_optimized_B3P86.xyz'
            # xyz_split(filename)


run()
