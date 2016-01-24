from molecule.molecule import Molecule
from openbabel import OBAtomAtomIter, OBConversion, OBMol
from molecule.atom import GetAtomicSymbol, Atom
from a_calculations.csm_calculations import exact_calculation

def xyz_split(filename):
    '''
    :param filename: name of xyz file containing multiple molecules, to be split
    :return: dictionary owith key: index of molecule and val: molecule string (params are not included)
    '''
    mol_dict={}
    index=0
    with open (filename, 'r') as file:
        for line in file:
            try:
                mol=""
                num= int(line)
                mol=line
                mol+=file.readline()
                for x in range(num):
                    mol+=file.readline()
                index+=1
                mol_dict[index]=mol
            except:
                 continue
    return mol_dict

def get_symmetries_list(filename):
    with open (filename, 'r') as file:
        args=[]
        for line in file:
            CSM=line.find("csm")+3
            IN=line.find("__IN__")
            arg=line[CSM:IN]
            cleaned=arg.replace(" ","")
            args.append(cleaned)
    return args





def run():
    filename=r'C:\Users\devora.witty\Sources\csm\Testing\testing for proteins\tests\mols_for_Itay\input\input_1_methane_csm\methane-test1.xyz'
    mol_dict= xyz_split(filename)
    filename=r'C:\Users\devora.witty\Sources\csm\Testing\testing for proteins\tests\mols_for_Itay\expected_output\expected_output_1_methane_csm\sym.txt'
    symmetries= get_symmetries_list(filename)
    for index in mol_dict:
        xyz=mol_dict[index]
        molecule=Molecule.read_string(xyz, "xyz")
#        molecule.find_equivalence_classes()
        for symmetry in symmetries:
            result=exact_calculation(symmetry, molecule)
            #print("initial coordinates:", molecule.atom_cords())
            print("symmetry:", symmetry)
            print ("outatoms:", result.outAtoms)
            print("csm:", result.csm)
            print("dir:", result.dir)
            print("perm:", result.perm)
            return



    #filename=r'C:\Users\devora.witty\Sources\csm\Testing\testing for proteins\tests\mols_for_Itay\input\input_4-clusters\W_Au12_optimized_B3P86.xyz'
    #xyz_split(filename)

run()