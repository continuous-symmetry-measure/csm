from csm.main.csm_run import run
import sys
from csm.molecule.normalizations import normalize_coords, de_normalize_coords
from argparse import ArgumentParser
from csm.input_output.arguments import _create_parser
from csm.molecule.molecule import Molecule
from csm.calculations.exact_calculations import exact_calculation
import numpy as np
from argparse import RawTextHelpFormatter
import logging
logger = logging.getLogger(__name__)

def get_normalization_type(args):
    parser = _create_parser()
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
    return normalizations, norm_file

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




def get_chain_perm(molecule, perm):
    '''
    finds the existing permutation between chains, in the result
    :return:
    '''
    chain_perm_dict={}
    for chain in molecule.chains:
        index= molecule.chains[chain][0]
        permuted_index= perm[index]
        for chain2 in molecule.chains:
            if permuted_index in molecule.chains[chain2]:
                chain_perm_dict[chain]=chain2
                break
    chain_perm=[]
    for chain in molecule.chains:
        permuted_index= chain_perm_dict[chain]
        chain_perm.append(permuted_index)

    return chain_perm

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
    molecule=result.molecule
    original_norm = molecule.norm_factor ** 2
    #coords = result.normalized_molecule_coords
    #symm = result.normalized_symmetric_structure
    coords=[atom.pos for atom in result.molecule.atoms]
    symm=result.symmetric_structure


    if norm_type == '0':  #standard
        return original_norm, original_csm
    if norm_type == '1':    #1 center of masses of the fragments
        coords = result.normalized_molecule_coords
        fragment_centers = get_fragment_centers(molecule.chains, coords, file)
        #norm = sum of distance (between center of mass and atom) squared
        norm=get_norm_by_distance_from_centers(coords, molecule.chains, fragment_centers)
        #divide by norm
        print_numdenom(file, original_csm, norm)
        return norm*original_norm, original_csm/norm
    if norm_type == '2':    #2 normalization according to symmetry of fragments, with existing perm
        fragment_centers = get_fragment_centers(molecule.chains, coords, file)
        #create a dummy molecule made up of atoms located at center of each mass
        coordinates_array=[fragment_centers[chain] for chain in molecule.chains]
        dummy = Molecule.molecule_from_coords(coordinates_array, molecule.chain_equivalences)
        #get chain permutation
        perm=get_chain_perm(molecule, result.perm)
        #run CSM using the perm
        new_result=exact_calculation(result.op_type, result.op_order, dummy, perm=perm)
        write_new_molecule(file, new_result)
        new_symm=result.symmetric_structure
        #(save s0, print the received CSM and the symmetric structure (ie of the mass centers) and the dir)
        #find normalization factor based on the above step
        coords=result.normalized_molecule_coords
        coordinates_dict={chain:new_symm[i] for i,chain in enumerate(molecule.chains)}
        norm = get_norm_by_distance_from_centers(coords, molecule.chains, coordinates_dict)
        print_numdenom(file, original_csm*original_norm, norm)
        return norm, original_csm*original_norm/norm

    if norm_type == '3':    #3 normalization according to symmetry of fragments, withOUT existing perm
        fragment_centers = get_fragment_centers(molecule.chains, coords, file)
        #create a dummy molecule made up of atoms located at center of each mass
        coordinates_array=[fragment_centers[chain] for chain in molecule.chains]
        dummy = Molecule.molecule_from_coords(coordinates_array, molecule.chain_equivalences)
        #run CSM
        new_result=exact_calculation(result.op_type, result.op_order, dummy, no_constraint=True, suppress_print=True)
        write_new_molecule(file, new_result)
        new_symm=result.symmetric_structure
        #(save s0, print the received CSM and the symmetric structure (ie of the mass centers) and the dir)
        #find normalization factor based on the above step
        coordinates_dict={chain:new_symm[i] for i,chain in enumerate(molecule.chains)}
        norm = get_norm_by_distance_from_centers(coords, molecule.chains, coordinates_dict)
        print_numdenom(file, original_csm*original_norm, norm)
        return norm, original_csm*original_norm/norm

    if norm_type == '4':    #4 normalization according to averages of approximation to symmetry of fragments
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
        numerator = 0
        for i in range(len(coords)):
            numerator += np.linalg.norm(coords[i] - symm[i])
        norm = len(coords)
        print_numdenom(file, numerator*100, norm)
        return norm, numerator*100 / norm


    if norm_type == '6': #6 Linear normalization
        #similar to standard csm but no squaring in numerator/denominator
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


def normrun(args=[]):
    if not args:
        args = sys.argv[1:]

    norm_types, norm_file = get_normalization_type(args)
    args=[x for x in args if x not in norm_types]  # remove the normalization argument

    result = run(args)


    if not set(norm_types).isdisjoint(('1','2','3','4')):
        if len(result.molecule.chains) <= 1:
            raise ValueError("Normalization types 1,2,3,4 are based on the molecule's fragments, "
                           "and the input molecule does not have multiple fragments.")

    file=None
    if norm_file:
        file = open(norm_file, 'w')
        write_coords(file, result.molecule.atoms, "Normalized coords", result.normalized_molecule_coords)

    for norm_type in norm_types:
        if file:
            file.write("\n***Normalization "+ norm_type +":***\n")
        print("--------")
        try:
            norm_factor, final_csm = normalize_csm(norm_type, result, file)
            print("Csm normalized with", normalization_dict[norm_type], "("+ norm_type+ ")", "is:", final_csm)
            print("Normalization factor is:", norm_factor)
        except Exception as e:
            print("FAILED to normalize csm with",  normalization_dict[norm_type], "("+ norm_type+ ")")
            print("Cause:", str(e))

    if file:
        file.close()


def run_no_return(args=[]):
    normrun(args)


if __name__ == '__main__':
    normrun(args=sys.argv[1:])
