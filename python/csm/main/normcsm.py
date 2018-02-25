from csm.calculations import ExactCalculation
from csm.calculations.data_classes import Operation
from csm.main.csm_run import run_calculation as csmrun
import sys
from csm.molecule.normalizations import normalize_coords, de_normalize_coords
from argparse import ArgumentParser
from csm.molecule.molecule import Molecule, MoleculeFactory
import numpy as np
from argparse import RawTextHelpFormatter
import logging
logger = logging.getLogger(__name__)


def exact_calculation(op_type, op_order, molecule, sn_max=8, keep_structure=False, perm=None, no_constraint=False, suppress_print=False, timeout=300, *args, **kwargs):
    ec= ExactCalculation(Operation.placeholder(op_type, op_order, sn_max), molecule, keep_structure, perm, no_constraint, timeout)
    ec.calculate()
    return ec.result


def old_create_parser():
    parser = ArgumentParser(usage="\ncsm type input_molecule output_file [additional arguments]", allow_abbrev=False)

    parser.add_argument('type',
                        # choices=c_symmetries + s_symmetries + ['cs', 'ci', 'ch'],
                        help='The type of operation: cs, ci, ch, cN, sN')
    parser.add_argument('input', help='Input molecule file')
    parser.add_argument('output', default='output.txt', help='Output file')

    # Optional arguments (their names start with --)

    # types of calculations (default is exact):
    calculation_type = parser.add_argument_group('Calculation Type (default is exact)')
    calculation_type.add_argument('--approx', action='store_true', default=False,
                                  help='use the approximate algorithm to estimate the CSM')
    calculation_type.add_argument('--trivial', action='store_true', default=False,
                                  help='CSM of identity perm, or, if chains, CSM of chain permutation with no atom permutation')


    parser.add_argument('--timeout', default=300,
                        help="Specify a timeout for CSM in seconds. Default is 5 minutes (300)", type=int)
    parser.add_argument('--sn-max', type=int, default=8, help='The maximal sn to try, relevant only for chirality')

    # general input/calculation arguments:
    # parser.add_argument('--ignore-hy', action='store_true', default=False, help='Ignore Hydrogen atoms in computations')
    input_type = parser.add_argument_group('Input Arguments')
    input_type.add_argument('--remove-hy', action='store_true', default=False,
                            help='Remove Hydrogen atoms, rebuild molecule without them, and compute')
    input_type.add_argument('--ignore-sym', action='store_true', default=False,
                            help='Ignore all atomic symbols, performing a purely geometric operation')
    input_type.add_argument('--use-mass', action='store_true', default=False,
                            help='Use the atomic masses to define center of mass')
    input_type.add_argument('--babel-bond', action='store_true', default=False, help='Let OpenBabel compute bonding')
    # parser.add_argument('--no-babel',  action='store_true', default=False, help='force suppress automatically using OpenBabel to compute bonds')
    input_type.add_argument('--use-sequence', action='store_true', default=False,
                            help='create equivalence class for pdb file using sequence information. Can\'t be used with --use-chains')
    input_type.add_argument('--use-chains', action='store_true', default=False,
                            help='When a molecule has chains, use them (affects trivial, approx)')
    input_type.add_argument('--read-fragments', action='store_true', default=False,
                            help='Read fragments from .mol or .pdb file as chains')


    # calculation arguments that only apply to exact:
    exact_args = parser.add_argument_group('Arguments for Exact Algorithm')
    exact_args.add_argument('--use-perm', type=str,
                            help='EXACT ONLY: Compute exact CSM for a single permutation')
    exact_args.add_argument('--keep-structure', action='store_true', default=False,
                            help='EXACT ONLY: Maintain molecule structure from being distorted in the exact calculation')
    exact_args.add_argument('--no-constraint', action='store_true', default=False,
                            help='EXACT ONLY: Do not use the constraints algorithm to traverse the permutation tree')

    # calculation arguments that only apply to approx
    # parser.add_argument('--use-dir', type=str,
    #                    help='Run the approx algorithm using predefined axes as the starting point')
    approx_args = parser.add_argument_group('Arguments for Approx Algorithm')
    approx_args.add_argument('--detect-outliers', action='store_true', default=False,
                             help="APPROX ONLY:Use outlier detection to improve guesses for initial directions in approx algorithm")
    approx_args.add_argument('--no-orthogonal', action='store_true', default=False,
                             help="APPROX ONLY:Don't add orthogonal directions to calculated directions")
    approx_args.add_argument('--fibonacci', type=int,
                             help="APPROX ONLY: Use fibonacci sphere to generate 50 starting directions")
    approx_args.add_argument('--use-best-dir', action='store_true', default=False,
                             help='APPROX ONLY:Only use the best direction')
    approx_args.add_argument('--many-chains', action='store_true', default=False,
                             help='APPROX ONLY: Use the new chains algorithm for many chains. Will automatically apply use-chains')
    approx_args.add_argument('--greedy', action='store_true', default=False,
                             help='APPROX ONLY: use the old greedy approx algorithm (no hungarian)')
    approx_args.add_argument('--dir', nargs=3, type=float,
                             help='run approximate algorithm using a specific starting direction')
    approx_args.add_argument('--selective', type=int,
                        help='Do a single iteration on many directions (use with --fibonacci), and then a full set of iterations only on the best k (default 10)')
    approx_args.add_argument('--statistics', type=str,
                        help='Print statistics about each direction to a file')
    approx_args.add_argument('--polar', action='store_true', default=False,
                        help="Print polar coordinates instead of cartesian coordinates")

    # output formatting and printing options
    out_args = parser.add_argument_group("Output Arguments")
    out_args.add_argument('--format', help='Use a specific input/output format')
    out_args.add_argument('--json-output', action='store_true', default=False,
                        help='Print output in json format to a file')
    out_args.add_argument('--print-local', action='store_true', default=False,
                        help='Print the local CSM (csm for each atom) in the output file')
    out_args.add_argument('--output-perms', action='store', default=None,
                        help='Writes all enumerated permutations to file')
    out_args.add_argument('--output-branches', action='store_true', default=False,
                        help='Writes all backtracking branches to the console')
    out_args.add_argument('--print-approx', action='store_true', default=False,
                        help='add some printouts to approx')
    out_args.add_argument('--print-denorm', action='store_true', default=False,
                        help='when printing the original molecule, print the denormalized coordinates')

    # defunct: no longer applied in code
    # parser.add_argument('--no-limit', action='store_true', default=False, help='Allows running program while ignoring computational complexity')
    # parser.add_argument('--babel-test', action='store_true', default=False, help="Test if the molecule is legal or not")
    # parser.add_argument('--time-only', action='store_true', default=False, help="Only print the time and exit")
    # parser.add_argument('--write-openu', action='store_true', default=False,
    #                    help='Write output in open university format')


    return parser

def old_process_arguments(parse_res):
    """
    Divides the parsed arguments (from argparse) into three dictionaries
    Args:
        parse_res: Result of argparse

    Returns:
        dictionary_args - input arguments, calculation arguments and output arguments

    """

    dictionary_args = {}
    dictionary_args["parallel"]=False
    # the first three positional arguments

    op = Operation(parse_res.type)

    dictionary_args['operation'] = op
    dictionary_args['op_type'] = op.type
    dictionary_args['op_order'] = op.order
    dictionary_args['op_name'] = op.name

    dictionary_args['in_file_name'] = parse_res.input

    dictionary_args['out_file_name'] = parse_res.output

    # optional arguments:
    dictionary_args['timeout'] = parse_res.timeout

    # types of calculations:
    dictionary_args['command'] = 'exact'
    if parse_res.approx:
        dictionary_args['command'] = 'approx'
    if parse_res.trivial:
        dictionary_args['command'] = 'trivial'

    # general input/calculation arguments:
    # dictionary_args['ignore_hy'] = parse_res.ignore_hy
    dictionary_args['remove_hy'] = parse_res.remove_hy
    dictionary_args['ignore_symm'] = parse_res.ignore_sym
    dictionary_args['sn_max'] = parse_res.sn_max
    dictionary_args['use_mass'] = parse_res.use_mass
    dictionary_args['babel_bond'] = parse_res.babel_bond
    dictionary_args['use_sequence'] = parse_res.use_sequence

    # if parse_res.use_sequence and parse_res.keep_structure:
    #    raise ValueError("--keep-structure and --use-sequence are mutually exclusive")



    # use chains and fragments
    dictionary_args['use_chains'] = parse_res.use_chains
    dictionary_args['read_fragments'] = parse_res.read_fragments

    if not dictionary_args['use_chains'] and parse_res.read_fragments:
        dictionary_args['use_chains'] = True

    # calculation arguments for exact only:
    if parse_res.use_perm:
        dictionary_args['perm_file_name'] = parse_res.use_perm

    dictionary_args['keep_structure'] = parse_res.keep_structure

    dictionary_args['no_constraint'] = parse_res.no_constraint


    # calculation arguments for approx only:
    dictionary_args['approx_algorithm'] = 'hungarian'
    if parse_res.many_chains:
        dictionary_args['use_chains'] = True
        dictionary_args['approx_algorithm'] = 'many-chains'
    if parse_res.greedy:
        dictionary_args['approx_algorithm'] = 'greedy'

    dictionary_args['detect_outliers'] = parse_res.detect_outliers
    dictionary_args['get_orthogonal'] = not parse_res.no_orthogonal
    dictionary_args['use_best_dir'] = parse_res.use_best_dir

    if parse_res.fibonacci is not None:
        dictionary_args["fibonacci"] = True
        dictionary_args["num_dirs"] = parse_res.fibonacci

    if parse_res.selective is not None:
        dictionary_args["selective"] = True
        dictionary_args["num_selected"] = parse_res.selective

    dir = parse_res.dir
    if dir:
        dictionary_args['dirs'] = [dir]

    # output arguments:
    dictionary_args['json_output'] = parse_res.json_output
    dictionary_args['print_approx'] = parse_res.print_approx
    dictionary_args['print_perms'] = parse_res.output_perms
    dictionary_args['print_branches'] = parse_res.output_branches
    dictionary_args['print_denorm'] = parse_res.print_denorm
    dictionary_args['format'] = parse_res.format
    dictionary_args['useformat'] = dictionary_args['format'] is not None
    if not dictionary_args['format']:
        # get input file extension
        dictionary_args['format'] = parse_res.input.split(".")[-1]

    # dictionary_args['write_openu'] = parse_res.write_openu
    dictionary_args['print_local'] = dictionary_args['calc_local'] = parse_res.print_local

    dictionary_args['perms_csv_name'] = parse_res.output_perms
    dictionary_args['polar']=parse_res.polar
    dictionary_args['stat_file_name']=parse_res.statistics
    return dictionary_args




def get_normalization_type(args):
    parser = old_create_parser()
    parser.formatter_class=RawTextHelpFormatter
    parser.usage = "\nnorm_csm normalization type input_molecule output_file [additional arguments]"
    norm_argument= parser.add_argument('normalization', default='0',
                        help = 'Types of normalization available:\n'
                               '0: standard normalization, according to centers of mass (without scaling)\n'
                               '1: normalization according to the center of mass of each fragment\n'
                               '2: normalization according to an approximation of the symmetric structure of the centers '
                               'of mass of each fragment, based on the solution permutation\n'
                               '3: normalization according to an approximation of the symmetric structure of the centers '
                               'of mass of each fragment, without using the solution permutation\n'
                               '4: normalization according to averages of approximation to symmetry of fragments\n'
                               '5: normalization according to number of atoms\n'
                               '6: linear normalization',
                        choices=['0', '1', '2', '3', '4', '5', '6'],
                        nargs='+', metavar="normalization"
                        )
    parser._actions.pop()
    parser._actions.insert(1, norm_argument)
    parser.add_argument('--output-norm', action='store', default=None,
                        help='Write debug information from normalization factors to a file')
    parsed_args = parser.parse_args(args)
    normalizations = parsed_args.normalization
    norm_file=parsed_args.output_norm
    if norm_file: #remove the norm file arguments
        args.remove('--output-norm')
        args.remove(norm_file)

    #TODO: add check that by perm is only if keep-structure or use-chains is applied
    if not set(normalizations).isdisjoint(('1', '2', '3', '4')):
        if not parsed_args.use_chains:
            raise ValueError("You selected a normalization type (1,2,3, or 4) that expects fragments, but did not select --use-chains")
    #TODO: add check that anything using fragments is either a pdb with chains, or includes a fragment file
    dictionary_args=old_process_arguments(parsed_args)
    return normalizations, norm_file,dictionary_args

#    dictionary_args['norm_output_file']= parse_res.output_norm

def get_fragments():
    pass

def get_fragment_centers(chains, positions, file):
    fragment_centers={}
    for chain in chains:
        fragment_centers[chain]=np.array([0.0, 0.0, 0.0])
        for index in chains[chain]:
            fragment_centers[chain]+=positions[index]
        fragment_centers[chain]/=len(chains[chain])

    if file:
        file.write("\nfragment centers:\n")
        for chainkey in fragment_centers:
            file.write("chain " + str(chainkey))
            file.write(str(fragment_centers[chainkey]))
            file.write("\n")
    return fragment_centers

#def divide_by_chain_centers(chains, positions):
#    fragment_centers = get_fragment_centers(chains, positions)
#    norm = 0
#    for chain in chains:
#        for index in chains[chain]:
#            norm += np.linalg.norm(positions[index] - fragment_centers[chain])
            # (molecule.Q[index] - fragment_centers[chain]) * (molecule.Q[index] - fragment_centers[chain])
#    return norm

def get_norm_by_distance_from_centers(coords, fragments, centers):
    '''
    :param coords: a list of coordinates whose distance from center will be measured
    :param fragments: a dictionary, key:fragment, value: array of indiced within coords
    :param centers: a dictionary of centers, key:fragment, value: coordinate of center
    :return:
    '''
    norm=0
    for fragment in fragments:
        for index in fragments[fragment]:
            ci=coords[index]
            cf=centers[fragment]
            norm += np.square((np.linalg.norm(ci - cf)))
    return norm




def write_new_molecule(file, result):
    if file:
        file.write("\n")
        write_coords(file, result.molecule.atoms, "dummy molecule coordinates")
        file.write("\n")
        write_coords(file, result.molecule.atoms, "symmetric structure", result.symmetric_structure)
        file.write("\n")
        file.write("\ndummy molecule's equivalence classes\n")
        file.write(str(result.molecule.equivalence_classes))
        file.write("\ncsm of dummy molecule\n")
        file.write(str(result.csm))
        file.write("\nperm\n")
        file.write(str([i+1 for i in result.perm]))
        file.write("\n\n")

def write_coords(file, atoms, header="", positions=None):
    file.write(Molecule.xyz_string(atoms, positions, header))

def print_numdenom(file, numerator, denominator):
    if file:
        file.write("\nNumerator: ")
        file.write(str(numerator))
        file.write("\nDenominator: ")
        file.write(str(denominator))
        file.write("\n")

def normalize_csm(norm_type, result, file):
    '''
    :param norm_type: the type of normalization factor
    :param result: we run this after having run CSM, so this is the result we received from running CSM
    :return: norm_factor, csm
    '''
    original_csm = result.csm
    result.molecule.create_Q() #make sure Q is up to date
    molecule=result.molecule
    original_norm = molecule.norm_factor ** 2
    #coords = result.normalized_molecule_coords
    #symm = result.normalized_symmetric_structure



    if norm_type == '0':  #standard
        return original_norm, original_csm
    if norm_type == '1':    #1 center of masses of the fragments
        coords = result.normalized_molecule_coords
        fragment_centers = get_fragment_centers(molecule.chains, coords, file)
        norm=get_norm_by_distance_from_centers(coords, molecule.chains, fragment_centers)
        print_numdenom(file, original_csm, norm)
        return norm*original_norm, original_csm/norm
    if norm_type == '2':    #2 normalization according to symmetry of fragments, with existing perm
        coords = result.molecule.Q
        fragment_centers = get_fragment_centers(molecule.chains, coords, file)
        #create a dummy molecule made up of atoms located at center of each mass
        coordinates_array=[fragment_centers[chain] for chain in molecule.chains]
        dummy = MoleculeFactory.dummy_molecule_from_coords(coordinates_array, molecule.chain_equivalences)
        #get chain permutation
        perm=result.chain_perm
        #run CSM using the perm
        new_result=exact_calculation(result.op_type, result.op_order, dummy, perm=perm)
        write_new_molecule(file, new_result)
        new_symm=new_result.symmetric_structure
        #find normalization factor based on the above step
        coordinates_dict={chain:new_symm[i] for i,chain in enumerate(molecule.chains)}
        norm = get_norm_by_distance_from_centers(coords, molecule.chains, coordinates_dict)
        print_numdenom(file, original_csm*original_norm, norm)
        return norm, original_csm*original_norm/norm

    if norm_type == '3':    #3 normalization according to symmetry of fragments, withOUT existing perm
        coords = result.molecule.Q
        fragment_centers = get_fragment_centers(molecule.chains, coords, file)
        #create a dummy molecule made up of atoms located at center of each mass
        coordinates_array=[fragment_centers[chain] for chain in molecule.chains]
        dummy = MoleculeFactory.dummy_molecule_from_coords(coordinates_array, molecule.chain_equivalences)
        #run CSM
        new_result=exact_calculation(result.op_type, result.op_order, dummy, suppress_print=True)
        write_new_molecule(file, new_result)
        new_symm=new_result.symmetric_structure
        #(save s0, print the received CSM and the symmetric structure (ie of the mass centers) and the dir)
        #find normalization factor based on the above step
        coordinates_dict={chain:new_symm[i] for i,chain in enumerate(molecule.chains)}
        norm = get_norm_by_distance_from_centers(coords, molecule.chains, coordinates_dict)
        print_numdenom(file, original_csm*original_norm, norm)
        return norm, original_csm*original_norm/norm

    if norm_type == '4':    #4 normalization according to averages of approximation to symmetry of fragments
        coords = result.molecule.Q
        symm = result.symmetric_structure
        #find center of mass of each fragment in the symmetric structure
        fragment_centers= get_fragment_centers(molecule.chains, symm, file)
        norm =get_norm_by_distance_from_centers(coords, molecule.chains, fragment_centers)
        #divide by norm
        print_numdenom(file, original_csm*original_norm, norm)
        return norm, original_csm*original_norm/norm

    if norm_type == '5': #5 normalization by number of atoms
        #atom factor validity can be tested by:
        #  multiplying normalized_coords and normalized_symm by x
        #  and verifying that the returned csm is also mutiplied by x
        coords = result.molecule.Q
        symm = result.symmetric_structure
        numerator = 0
        for i in range(len(coords)):
            numerator += np.linalg.norm(coords[i] - symm[i])
        norm = len(coords)
        print_numdenom(file, numerator*100, norm)
        return norm, numerator*100 / norm


    if norm_type == '6': #6 Linear normalization
        #similar to standard csm but no squaring in numerator/denominator
        coords = result.molecule.Q
        symm = result.symmetric_structure
        numerator = norm = 0
        for i in range(len(coords)):
            numerator += np.linalg.norm(coords[i] - symm[i])
            norm += np.linalg.norm(coords[i])  # we assume that the center of mass of the whole molecule is (0,0,0).
        print_numdenom(file, numerator*100, norm)
        return norm, numerator*100 / norm



normalization_dict={
    "0": "standard normalization",
    "1": "fragment center",
    "2": "fragment center symmetric structure with perm",
    "3": "fragment center symmetric structure (no perm)",
    "4": "averages of approximation of symmetric centers",
    "5": "number of atoms",
    "6": "linear normalization"
}


def run(args=[]):

    if not args:
        args = sys.argv[1:]

    norm_types, norm_file, parsed_args = get_normalization_type(args)
    args=[x for x in args if x not in norm_types]  # remove the normalization argument

    result = csmrun(parsed_args)


    if not set(norm_types).isdisjoint(('1','2','3','4')):
        if len(result.molecule.chains) <= 1:
            raise ValueError("Normalization types 1,2,3,4 are based on the molecule's fragments, "
                           "and the input molecule does not have multiple fragments.")

    normalization_results={}

    file=None
    try:
        if norm_file:
            file = open(norm_file, 'w')
            write_coords(file, result.molecule.atoms, "Normalized coords", result.normalized_molecule_coords)

        for norm_type in norm_types:
            if file:
                file.write("\n***Normalization "+ norm_type +":***\n")
            print("--------")
            try:
                norm_factor, final_csm = normalize_csm(norm_type, result, file)
                normalization_results[norm_type]=(norm_factor, final_csm)
                print("Csm normalized with", normalization_dict[norm_type], "("+ norm_type+ ")", "is:", final_csm)
                print("Normalization factor is:", norm_factor)
            except IOError as e:
                print("FAILED to normalize csm with",  normalization_dict[norm_type], "("+ norm_type+ ")")
                print("Cause:", str(e))
        return normalization_results

    finally:
        if file:
            file.close()


def run_norm_no_return(args=[]):
    run(args)

if __name__ == '__main__':
    run(args=sys.argv[1:])



