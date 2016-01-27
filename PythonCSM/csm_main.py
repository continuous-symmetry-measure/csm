import sys
import logging
from input_output.arguments import get_split_arguments, get_arguments
from molecule.molecule import Molecule
from a_calculations.csm_calculations_data import CSMCalculationsData
from a_calculations.csm_calculations import approx_calculation, exact_calculation, local_calculation

MINDOUBLE = 1e-8
APPROX_RUN_PER_SEC = 8e4
sys.setrecursionlimit(10000)

logger = None

def init_logging(log_file_name=None):
    global logger

    if log_file_name:
        logging.basicConfig(filename=log_file_name, level=logging.DEBUG,
                            format='[%(asctime)-15s] [%(levelname)s] [%(name)s]: %(message)s')
    else:
        logging.basicConfig(level=logging.ERROR)
    logger = logging.getLogger("csm")

def run_csm(args={}, print_output=True):

    #initialize
    mol_args, calc_args, out_args=get_split_arguments(args)
    mol=Molecule.read(**mol_args)

    #backwards compatibility
    csm_args = get_arguments(args)
    csm_args['molecule']=mol
    #csm.SetCSMOptions(csm_args)
    cppdata = CSMCalculationsData(csm_args)
    calc_args['cppdata']=cppdata

    #logging:
    init_logging(**out_args)

    #run actual calculation
    #approx:
    if calc_args['approx']:
        approx_calculation(**calc_args)
    #local:

    #exact:




    #print results
    #opName, csm, scalingfactor=dmin, dir, equivalence classes, molecule, localCSM, chMinOrder, perm
    '''
            if csm_args['printLocal']:
            if csm_args['type'] == 'CH':
                data.opOrder = result.chMinOrder
            local_res = csm.ComputeLocalCSM(data)
            result.localCSM = local_res.localCSM

        process_results(result, csm_args)

        if print_output:
            print_all_output(result, csm_args)

        return result
    finally:
        try:
            csm_args['outFile'].close()
            csm_args['outPermFile'].close()
        except:
            pass
    '''


def read_dir_file(f):
    """
    Reads a symmetry direction file
    :param f: File object
    :return: (x,y,z) of the symmetry axis
    """
    line = f.readline().split()
    result = (float(line[0]), float(line[1]), float(line[2]))
    return result

def read_perm_file(f):
    """
    Reads a permutation
    :param f: File object
    :return: permutation as a list of numbers

    Check that the permutation is legal, raise ValueError if not
    """
    line = f.readline().split()
    used = []
    for i in range(len(line)):
        used.append(False)

    result = []
    for num_str in line:
        try:
            num = int(num_str)
        except ValueError:
            raise ValueError("Invalid permutation")
        if num < 1 or num > len(line) or used[num-1]:
            raise ValueError("Invalid permutation")
        result.append(num-1)
        used[num-1] = True
    return result



if __name__ == '__main__':
    results = run_csm(args=sys.argv[1:])