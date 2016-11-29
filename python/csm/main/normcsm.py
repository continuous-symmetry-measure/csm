from csm.main.csm_run import run
import sys
from csm.molecule.normalizations import normalize_coords, de_normalize_coords
from argparse import ArgumentParser
from csm.input_output.arguments import _create_parser
from csm.molecule.molecule import Molecule
from csm.calculations.exact_calculations import exact_calculation
import numpy as np

def get_normalization_type(args):
    parser = _create_parser()
    parser.usage = "\nnormalization csm type input_molecule output_file [additional arguments]"
    norm_argument= parser.add_argument('normalization', default='standard',
                        help = 'Types of normalizations: standard, fragment_center, fragment_perm, fragment_symm, symmetry_center, atom_number, linear_csm',
                        choices=['standard', 'atom_number', 'fragment_center', 'symmetry_center', 'fragment_symm', 'fragment_perm', 'linear_csm'],
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
    :param coords: a list of coordinates
    :param fragments: a dictionary of fragments, whose values are indexes within coordinated
    :param centers: a dictionary of centers, whose keys match the keys in fragments
    :return:
    '''
    norm=0
    for fragment in fragments:
        for index in fragments[fragment]:
            norm += np.square((np.linalg.norm(coords[index] - centers[fragment])))
    return norm


def atom_number_factor(coords, symm, original_norm):
    Pk_Qk = 0
    for i in range(len(coords)):
        Pk_Qk += np.linalg.norm(coords[i] - symm[i])
    norm_factor = len(coords) / original_norm
    new_csm = Pk_Qk / norm_factor
    return norm_factor, new_csm


def get_chain_perm(molecule, perm):
    chain_perm_dict={}
    for chain in molecule.chains:
        index= molecule.chains[chain][0]
        permuted_index= perm[index]
        for chain2 in molecule.chains:
            if permuted_index in molecule.chains[chain2]:
                chain_perm_dict[chain]=chain2
                break
    chain_perm=[]
    for i in range(len(molecule.reversechainkeys)):
        chain = molecule.reversechainkeys[i]
        permuted= chain_perm_dict[chain]
        permuted_index= molecule.chainkeys[permuted]
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
    original_norm = molecule.norm_factor
    normalized_coords = result.normalized_molecule_coords
    normalized_symm = result.normalized_symmetric_structure


    if norm_type == 'standard':  #0
        return original_norm, original_csm
    if norm_type == 'fragment_center':    #1 center of masses of the fragments
        #find the center of mass of each fragment
        fragment_centers= get_fragment_centers(molecule.chains, normalized_coords)
        #norm = sum of distance (between center of mass and atom) squared
        norm=get_norm_by_distance_from_centers(normalized_coords, molecule.chains, fragment_centers)
        #divide by norm
        return norm*original_norm, original_csm/norm
    if norm_type == 'fragment_perm':    #2 normalization according to symmetry of fragments, with existing perm
        #find center of mass of each fragment
        fragment_centers= get_fragment_centers(molecule.chains, normalized_coords)
        #create a dummy molecule made up of atoms located at center of each mass
        coordinates_array=[fragment_centers[molecule.reversechainkeys[i]] for i in range(len(molecule.reversechainkeys))]
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
        coordinates_dict={molecule.reversechainkeys[i]:new_symm[i] for i in range(len(molecule.reversechainkeys))}
        norm = get_norm_by_distance_from_centers(normalized_coords, molecule.chains, coordinates_dict)

        return norm * original_norm, original_csm / norm

    if norm_type == 'fragment_symm':    #3 normalization according to symmetry of fragments, withOUT existing perm
        #find center of mass of each fragment
        fragment_centers= get_fragment_centers(molecule.chains, normalized_coords)
        #create a dummy molecule made up of atoms located at center of each mass
        coordinates_array=[fragment_centers[molecule.reversechainkeys[i]] for i in range(len(molecule.reversechainkeys))]
        dummy = Molecule.molecule_from_coords(coordinates_array, molecule.chain_equivalences)
        #run CSM
        new_result=exact_calculation(result.op_type, result.op_order, dummy)
        #receive s0
        s0= new_result.csm
        new_symm=result.symmetric_structure
        #(save s0, print the received CSM and the symmetric structure (ie of the mass centers) and the dir)
        #find normalization factor based on the above step
        coordinates_dict={molecule.reversechainkeys[i]:new_symm[i] for i in range(len(molecule.reversechainkeys))}
        norm = get_norm_by_distance_from_centers(normalized_coords, molecule.chains, coordinates_dict)

        return norm * original_norm, original_csm / norm

    if norm_type == 'symmetry_center':    #4 normalization according to averages of approximation to symmetry of fragments
        #find center of mass of each fragment in the symmetric structure
        fragment_centers= get_fragment_centers(molecule.chains, normalized_symm)
        norm=get_norm_by_distance_from_centers(normalized_coords, molecule.chains, fragment_centers)
        #divide by norm
        return norm*original_norm, original_csm/norm

    if norm_type == 'atom_number': #5 normalization by number of atoms
        #note-- atom_number_factor was refactored as a separate equation because
        #it is possible to test its validity by sending 5*coord, 5*symm and seeing that
        #indeed the CSM increases times 5
        return atom_number_factor(normalized_coords, normalized_symm, original_norm)

    if norm_type == 'linear_csm':
        numerator = denominator = 0
        for i in range(len(normalized_coords)):
            numerator += np.linalg.norm(normalized_coords[i] - normalized_symm[i])
            denominator += np.linalg.norm(normalized_coords[i])  # we assume that the center of mass of the whole molecule is (0,0,0).
        return denominator * original_norm, numerator / denominator






def normrun(args=[]):
    if not args:
        args = sys.argv[1:]

    norm_types = get_normalization_type(args)
    args=[x for x in args if x not in norm_types]  # remove the normalization argument

    result = run(args)

    for norm_type in norm_types:
        print("--------")
        norm_factor, final_csm = normalize_csm(norm_type, result)
        print("Csm normalized with", norm_type, "method is:", final_csm)
        print("Normalization factor is:", norm_factor)


def run_no_return(args=[]):
    normrun(args)


if __name__ == '__main__':
    normrun(args=sys.argv[1:])
