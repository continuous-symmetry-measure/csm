from csm.main.csm_run import run
import sys
from csm.molecule.normalizations import normalize_coords, de_normalize_coords
from argparse import ArgumentParser
from csm.input_output.arguments import _create_parser
import numpy as np

def get_normalization_type(args):
    parser = _create_parser()
    parser.usage = "\ncsm type input_molecule output_file normalization [additional arguments]"
    parser.add_argument('normalization', default='standard', help='The type of normalization to apply',
                        choices=['standard', 'atom_number', 'fragment_mass_center', 'symmetric_fragment_mass_center'],
                        )
    parsed_args = parser.parse_args(args)
    normalization = parsed_args.normalization
    return normalization

def get_fragment_centers(chains, positions):
    fragment_centers={}
    for chain in chains:
        fragment_centers[chain]=np.array([0.0, 0.0, 0.0])
        for index in chains[chain]:
            fragment_centers[chain]+=positions[index]
        fragment_centers[chain]/=len(chains[chain])
    return fragment_centers

def divide_by_chain_centers(chains, positions):
    fragment_masses = get_fragment_centers(chains, positions)
    norm = 0
    for chain in chains:
        for index in chains[chain]:
            norm += np.linalg.norm(positions[index] - fragment_masses[chain])
            # (molecule.Q[index] - fragment_masses[chain]) * (molecule.Q[index] - fragment_masses[chain])
    return norm


def normalize_csm(norm_type, result):
    '''
    :param norm_type: the type of normalization factor
    :param result: we run this after having run CSM, so this is the result we received from running CSM
    :return:
    '''
    original_csm = result.csm
    molecule=result.molecule
    original_norm = molecule.norm_factor #this is the original normalization factor
    denormed_csm = original_csm * original_norm #the result was divided by that factor, so we undo that by multiplication

    if norm_type == 'standard':
        return original_norm, original_csm
    if norm_type == 'atom_number':
        return len(molecule), denormed_csm / len(molecule)
    if norm_type == 'fragment_mass_center':
        norm=divide_by_chain_centers(molecule.chains, molecule.Q)
        return norm, denormed_csm/norm
    if norm_type == 'symmetric_fragment_mass_center':
        norm=divide_by_chain_centers(molecule.chains, result.symmetric_structure)
        return norm, denormed_csm/norm



def normrun(args=[]):
    if not args:
        args = sys.argv[1:]

    norm_type = get_normalization_type(args)
    args.pop(3)  # remove the normalization argument

    result = run(args)

    norm_factor, final_csm = normalize_csm(norm_type, result)
    print("normalization factor is:", norm_factor)
    print("Csm normalized with", norm_type, "method is:", final_csm)


def run_no_return(args=[]):
    normrun(args)


if __name__ == '__main__':
    normrun(args=sys.argv[1:])
