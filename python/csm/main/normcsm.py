from csm.main.csm_run import run
import sys
from csm.molecule.normalizations import normalize_coords, de_normalize_coords
from argparse import ArgumentParser
from csm.input_output.arguments import _create_parser
from csm.molecule.molecule import Molecule
from csm.calculations.exact_calculations import exact_calculation
import numpy as np
from argparse import RawTextHelpFormatter


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
    parsed_args = parser.parse_args(args)
    normalization = parsed_args.normalization

    #TODO: add check that by perm is only if keep-structure or use-chains is applied
    #TODO: add check that anything using fragments is either a pdb with chains, or includes a fragment file
    return normalization

def get_fragments():
    pass

def get_fragment_centers(chains, positions):
    fragment_centers={}
    for chain in chains:
        fragment_centers[chain]=np.array([0.0, 0.0, 0.0])
        for index in chains[chain]:
            fragment_centers[chain]+=positions[index]
        fragment_centers[chain]/=len(chains[chain])
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
            norm += np.square((np.linalg.norm(coords[index] - centers[fragment])))
    return norm


def atom_number_factor(coords, symm):
    Pk_Qk = 0
    for i in range(len(coords)):
        Pk_Qk += np.linalg.norm(coords[i] - symm[i])
    norm_factor = len(coords)
    new_csm = Pk_Qk / norm_factor
    return norm_factor, new_csm


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
        permuted= chain_perm_dict[chain]
        permuted_index= molecule.chains.index_map[permuted]
        chain_perm.append(permuted_index)

    return chain_perm

def normalize_csm(norm_type, result):
    '''
    :param norm_type: the type of normalization factor
    :param result: we run this after having run CSM, so this is the result we received from running CSM
    :return:
    '''
    original_csm = result.csm
    molecule=result.molecule
    original_norm = molecule.norm_factor ** 2
    normalized_coords = result.normalized_molecule_coords
    normalized_symm = result.normalized_symmetric_structure


    if norm_type == '0':  #standard
        return original_norm, original_csm
    if norm_type == '1':    #1 center of masses of the fragments
        #find the center of mass of each fragment
        fragment_centers= get_fragment_centers(molecule.chains, normalized_coords)
        #norm = sum of distance (between center of mass and atom) squared
        norm=get_norm_by_distance_from_centers(normalized_coords, molecule.chains, fragment_centers)
        #divide by norm
        return norm*original_norm, original_csm/norm
    if norm_type == '2':    #2 normalization according to symmetry of fragments, with existing perm
        #find center of mass of each fragment
        fragment_centers= get_fragment_centers(molecule.chains, normalized_coords)
        #create a dummy molecule made up of atoms located at center of each mass
        coordinates_array=[fragment_centers[chain] for chain in molecule.chains]
        dummy = Molecule.molecule_from_coords(coordinates_array, molecule.chain_equivalences)
        #get chain permutation
        perm=get_chain_perm(molecule, result.perm)
        #run CSM using the perm
        new_result=exact_calculation(result.op_type, result.op_order, dummy, perm=perm)
        #receive s0
        s0= new_result.csm
        new_symm=result.symmetric_structure
        #(save s0, print the received CSM and the symmetric structure (ie of the mass centers) and the dir)
        #find normalization factor based on the above step

        coordinates_dict={chain:new_symm[i] for i,chain in enumerate(molecule.chains)}
        norm = get_norm_by_distance_from_centers(normalized_coords, molecule.chains, coordinates_dict)

        return norm * original_norm, original_csm / norm

    if norm_type == '3':    #3 normalization according to symmetry of fragments, withOUT existing perm
        #find center of mass of each fragment
        fragment_centers= get_fragment_centers(molecule.chains, normalized_coords)
        #create a dummy molecule made up of atoms located at center of each mass
        coordinates_array=[fragment_centers[chain] for chain in molecule.chains]
        dummy = Molecule.molecule_from_coords(coordinates_array, molecule.chain_equivalences)
        #run CSM
        new_result=exact_calculation(result.op_type, result.op_order, dummy, no_constraint=True, suppress_print=True)
        #receive s0
        s0= new_result.csm
        new_symm=result.symmetric_structure
        #(save s0, print the received CSM and the symmetric structure (ie of the mass centers) and the dir)
        #find normalization factor based on the above step
        coordinates_dict={chain:new_symm[i] for i,chain in enumerate(molecule.chains)}
        norm = get_norm_by_distance_from_centers(normalized_coords, molecule.chains, coordinates_dict)

        return norm * original_norm, original_csm / norm

    if norm_type == '4':    #4 normalization according to averages of approximation to symmetry of fragments
        #find center of mass of each fragment in the symmetric structure
        fragment_centers= get_fragment_centers(molecule.chains, normalized_symm)
        norm=get_norm_by_distance_from_centers(normalized_coords, molecule.chains, fragment_centers)
        #divide by norm
        return norm*original_norm, original_csm/norm

    if norm_type == '5': #5 normalization by number of atoms
        #note-- atom_number_factor was refactored as a separate equation because
        #it is possible to test its validity by sending 5*coord, 5*symm and seeing that
        #indeed the CSM increases times 5
        return atom_number_factor(normalized_coords, normalized_symm)

    if norm_type == '6': #6 Linear normalization
        #similar to standard csm but no squaring in numerator/denominator
        numerator = denominator = 0
        for i in range(len(normalized_coords)):
            numerator += np.linalg.norm(normalized_coords[i] - normalized_symm[i])
            denominator += np.linalg.norm(normalized_coords[i])  # we assume that the center of mass of the whole molecule is (0,0,0).
        return denominator * np.sqrt(original_norm), numerator / denominator



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

    norm_types = get_normalization_type(args)
    args=[x for x in args if x not in norm_types]  # remove the normalization argument

    result = run(args)

    for norm_type in norm_types:
        print("--------")
        try:
            norm_factor, final_csm = normalize_csm(norm_type, result)
            print("Csm normalized with", normalization_dict[norm_type], "("+ norm_type+ ")", "is:", final_csm)
            print("Normalization factor is:", norm_factor)
        except:
            print("FAILED to normalize csm with",  normalization_dict[norm_type], "("+ norm_type+ ")")


def run_no_return(args=[]):
    normrun(args)


if __name__ == '__main__':
    normrun(args=sys.argv[1:])
