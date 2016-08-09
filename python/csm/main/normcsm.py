from csm.main.csm_run import run
import sys
from csm.molecule.normalizations import normalize_coords, de_normalize_coords
from argparse import ArgumentParser
from csm.input_output.arguments import _create_parser

def get_normalization_type(args):
    parser=_create_parser()
    parser.usage="\ncsm type input_molecule output_file normalization [additional arguments]"
    parser.add_argument('normalization', default='standard', help='The type of normalization to apply',
                        choices=['standard', 'david'],
                        )
    parsed_args = parser.parse_args(args)
    normalization= parsed_args.normalization
    return normalization

def david_norm(csm, molecule):
    return csm/len(molecule)

def normalize_csm(csm, norm_type, result):
    if norm_type=='david':
        return david_norm(csm, result.molecule)


def normrun(args=[]):
    if not args:
        args = sys.argv[1:]

    norm_type=get_normalization_type(args)
    args.pop(3) #remove the normalization argument

    result= run(args)
    molecule =result.molecule
    original_norm = molecule.norm_factor
    original_csm= result.csm
    denormed_csm= original_csm / original_norm
    final_csm= normalize_csm(denormed_csm, norm_type, result)
    print("Csm normalized with", norm_type, "method is:", final_csm)

    #modifications made to csm in original csm code:
    #(fabs(100 * (1.0 - csm / op_order)) (in fast_calc)



def normrun_no_return(args=[]):
    normrun(args)

if __name__ == '__main__':
    normrun(args=sys.argv[1:])

