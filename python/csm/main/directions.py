import sys
import os
import random
import numpy as np
from argparse import RawTextHelpFormatter

from csm.calculations.approx.dirs import dirs_orthogonal
from csm.calculations.constants import MAXDOUBLE
from csm.calculations.exact_calculations import exact_calculation, CSMValueError
from csm.calculations.approx.main import approx_calculation
from csm.input_output.arguments import _create_parser, get_split_arguments
from csm.molecule.molecule import Molecule, MoleculeFactory

__author__ = 'Devora Witty'

'''
accessory program for running different kinds of direction finding algorithms with the CSM approx algorithm
'''


def direction_parser():
    parser = _create_parser()
    parser.formatter_class = RawTextHelpFormatter
    parser.usage = "\ndirection_test direction_choice type input_molecule output_file [additional arguments]"
    dir_argument = parser.add_argument('direction_choice',
                                       help='Types of direction choice available:\n'
                                            'user-input: user inputs the dirs to be chosen with "--dirs file" \n'
                                            'exact-structure: program runs exact with keep structure and uses that direction for approx\n'
                                            'greedy-first: program runs standard greedy approx and uses final result direction for hungarian\n'
                                            'random-k: K directions are chosen randomly. Default is 10, user can specify --k K.\n'
                                            'cube-corners: directions to corners of cube whose center is mass center\n'
                                            'atom-vectors: directions are vectors from center to each atom\n'
                                            'atom-vector-orth: like atom-vectors, but also includes orthogonal directions to those chosen \n'
                                       ,
                                       choices=['user-input', 'exact-structure', 'greedy-first', 'random-k',
                                                'cube-corners', 'atom-vectors', 'atom-vectors-orth'],
                                       metavar="direction_choice"
                                       )
    parser._actions.pop()
    parser._actions.insert(1, dir_argument)
    parser.add_argument('--k', type=int, default=10, help='Number of random directions to try for random-k')
    parser.add_argument('--dirs-file', type=str,
                        help='File address of file with list of dirs for use-input')
    return parser


def choose_directions(direction_choice, dirs_file, k, molecule, csm_args):
    dirs = []
    op_type=csm_args['op_type']
    op_order= csm_args['op_order']
    if direction_choice == 'user-input':
        if not os.path.isfile(dirs_file):
            raise ValueError(
                "If direction type is user-input you must specify a valid dirs file using the --dirs-file flag")
        with open(dirs_file, 'r') as file:
            for line in file:
                dir = [float(x) for x in line.split()]
                if len(dir) != 3:
                    raise ValueError("Directions file is not correctly formatted")
                dirs.append(dir)

    if direction_choice == 'exact-structure':
        result=exact_calculation(op_type, op_order, molecule, keep_structure=True, suppress_print=True)
        print("Running exact calculation with keep-structure yielded CSM", result.csm)
        dirs.append(result.dir)

    if direction_choice == 'greedy-first':
        result=approx_calculation(op_type, op_order, molecule, approx_algorithm='greedy')
        print("Running approx calculation with greedy algorithm yielded CSM", result.csm)
        dirs.append(result.dir)

    if direction_choice == 'random-k':
        seed = random.randrange(sys.maxsize)
        rng = random.Random(seed)
        print("Seed was:", seed)

        for i in range(k):
            dir = []
            for j in range(3):
                dir.append(rng.uniform(-1, 1))
            dirs.append(dir)

    if direction_choice == 'cube-corners':
        pass

    if direction_choice == 'atom-vectors' or 'atom-vectors-orth':
        for atom in molecule.atoms:
            dir = [atom.pos[0] - molecule.center_of_mass[0],
                  atom.pos[1] - molecule.center_of_mass[1],
                  atom.pos[2] - molecule.center_of_mass[2]]
            dir /= np.linalg.norm(dir)
            dirs.append(dir)

    if direction_choice == 'atom-vectors-orth':
        dirs=dirs_orthogonal(dirs)

    return dirs


def run(args=[]):
    if not args:
        args = sys.argv[1:]
    if '--approx' not in args:
        args.append("--approx")
    if "--dir" in args:
        raise ValueError("The argument --dir cannot be used in conjunction with test_directions")

    parser = direction_parser()
    parsed_args = parser.parse_args(args)

    direction_choice = parsed_args.direction_choice
    dirs_file = parsed_args.dirs_file
    k = parsed_args.k

    del args[0]
    csm_args = get_split_arguments(args)
    # csm_args['op_type'], csm_args['op_order'], molecule, dir=dir,

    molecule  = MoleculeFactory.from_file(**csm_args)
    dirs = choose_directions(direction_choice, dirs_file, k, molecule, csm_args)

    best_csm=MAXDOUBLE
    best_dir=None

    for dir in dirs:
        csm_args['molecule'] = molecule.copy()
        csm_args['dirs'] = [dir]
        if True:
            try:
                result = approx_calculation(**csm_args)
                print("For initial direction", dir, ", the CSM value found was", result.csm, "with a final dir of",
                  result.dir)
                if result.csm<best_csm:
                    best_csm=result.csm
                    best_dir=dir
            except CSMValueError as e:
                result=e.CSMState
                print("***ERROR*** For initial direction", dir, ", the CSM value found was", result.csm, "with a final dir of",
                      result.dir)

    print("The best csm:", best_csm, "was found with the initial dir", best_dir)


if __name__ == '__main__':
    run(args=sys.argv[1:])
