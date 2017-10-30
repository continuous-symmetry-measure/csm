import sys
import os
import random
import numpy as np
from math import radians, cos, sin, pi
from argparse import RawTextHelpFormatter

from csm import __version__
from csm.calculations.approx.dirs import dirs_orthogonal
from csm.calculations.basic_calculations import check_perm_cycles, CSMState
from csm.calculations.constants import MAXDOUBLE
from csm.calculations.exact_calculations import exact_calculation, CSMValueError
from csm.calculations.approx.main import approx_calculation
from csm.input_output.arguments import _create_parser, get_split_arguments
from csm.input_output.formatters import format_CSM
from csm.input_output.writers import print_results
from csm.molecule.molecule import Molecule, MoleculeFactory, MoleculeReader

__author__ = 'Devora Witty'

'''
accessory program for running different kinds of direction finding algorithms with the CSM approx algorithm
'''

choice_dict = {
    '0': 'user-input',
    '1': 'exact-structure',
    '2': 'greedy-first',
    '3': 'random-k',
    '4': 'cube-corners',
    '5': 'atom-vectors',
    '6': 'atom-vectors-orth'
}


def direction_parser():
    parser = _create_parser()
    parser.formatter_class = RawTextHelpFormatter
    parser.usage = "\ndirection_test direction_choice dir_output type input_molecule output_file [additional arguments]"

    dir_argument = parser.add_argument('dir_output', default='dir-output.txt', help='Output file for directions')
    parser._actions.pop()
    parser._actions.insert(1, dir_argument)

    dir_argument = parser.add_argument('direction_choice',
                                       help='Types of direction choice available:\n'
                                            '0: user-input- user inputs the dirs to be chosen with "--dirs-file" \n'
                                            '1: exact-structure- program runs exact with keep structure and uses that direction for approx\n'
                                            '2: greedy-first- program runs standard greedy approx and uses final result direction for hungarian\n'
                                            '3: random-k- K directions are chosen randomly. Default is 10, user can specify --k K.\n'
                                            '4: cube-corners- directions to corners of cube whose center is mass center\n'
                                            '5: atom-vectors- directions are vectors from center to each atom\n'
                                            '6: atom-vector-orth- like atom-vectors, but also includes orthogonal directions to those chosen \n'
                                       ,
                                       choices=['0', '1', '2', '3', '4', '5', '6'],
                                       # nargs='+',
                                       metavar="direction_choice"
                                       )
    parser._actions.pop()
    parser._actions.insert(1, dir_argument)

    parser.add_argument('--num-rand', type=int, default=10, help='Number of random directions to try for random-k')
    parser.add_argument('--dirs-file', type=str,
                        help='File address of file with list of dirs for use-input')
    parser.add_argument('--seed', type=str,
                        help='If you\'d like to reproduce a run of random-k with a given seed')
    return parser


class PrintClass:
    file = ""

    @staticmethod
    def set(filename):
        PrintClass.file = filename

    @staticmethod
    def my_print(*strings, print_flag=False):
        if print_flag:
            print(*strings)
        with open(PrintClass.file, 'a') as file:
            for string in strings:
                file.write(str(string))
                file.write(" ")
            file.write("\n")


def choose_directions(direction_choice, molecule, csm_args, dirs_file, k, seed):
    dirs = []
    op_type = csm_args['op_type']
    op_order = csm_args['op_order']
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
        result = exact_calculation(op_type, op_order, molecule.copy(), keep_structure=True, suppress_print=True)
        PrintClass.my_print("\tRunning exact calculation with keep-structure yielded CSM", format_CSM(result.csm))
        dirs.append(result.dir)

    if direction_choice == 'greedy-first':
        result = approx_calculation(op_type, op_order, molecule.copy(), approx_algorithm='greedy')
        PrintClass.my_print("\tRunning approx calculation with greedy algorithm yielded CSM", format_CSM(result.csm))
        dirs.append(result.dir)

    if direction_choice == 'random-k':
        if not seed:
            seed = random.randrange(sys.maxsize)
        rng = random.Random(seed)
        PrintClass.my_print("\tSeed was:", seed)

        for i in range(k):
            dir = []
            for j in range(3):
                dir.append(rng.uniform(-1, 1))
            dirs.append(dir)

    if direction_choice == 'cube-corners':
        center = molecule.center_of_mass
        corner = np.add(center, [.5, .5, .5])
        dirs.append(corner)
        corner = np.add(center, [.5, .5, -.5])
        dirs.append(corner)
        corner = np.add(center, [.5, -.5, .5])
        dirs.append(corner)
        corner = np.add(center, [.5, -.5, -.5])
        dirs.append(corner)
        corner = np.add(center, [-.5, .5, .5])
        dirs.append(corner)
        corner = np.add(center, [-.5, .5, -.5])
        dirs.append(corner)
        corner = np.add(center, [-.5, -.5, .5])
        dirs.append(corner)
        corner = np.add(center, [-.5, -.5, -.5])
        dirs.append(corner)

    if direction_choice in ['atom-vectors', 'atom-vectors-orth']:
        for atom in molecule.atoms:
            dir = [atom.pos[0] - molecule.center_of_mass[0],
                   atom.pos[1] - molecule.center_of_mass[1],
                   atom.pos[2] - molecule.center_of_mass[2]]
            dir /= np.linalg.norm(dir)
            dirs.append(dir)

    if direction_choice == 'atom-vectors-orth':

        dirs = dirs_orthogonal(dirs)

    return dirs


def run_dir(dir, csm_args, molecule):
    csm_args['molecule'] = molecule.copy()
    csm_args['dirs'] = [dir]
    try:
        result = approx_calculation(**csm_args)
        PrintClass.my_print("\tFor initial direction", dir,
                            " the CSM value found was %s", format_CSM(result.csm),
                            " with a final dir of", result.dir)
    except CSMValueError as e:
        result = e.CSMState
        PrintClass.my_print("\t***FAILED TO FIND CSM*** For initial direction", dir,
                            " the CSM value found was", format_CSM(result.csm),
                            "with a final dir of", result.dir)

    falsecount, num_invalid, cycle_counts = check_perm_cycles(result.perm, result.op_order,
                                                              result.op_type)

    if False:
        PrintClass.my_print(
            "\tThe permutation found contains %d invalid %s. %.2lf%% of the molecule's atoms are in legal cycles" % (
                falsecount, "cycle" if falsecount == 1 else "cycles",
                100 * (len(result.molecule) - num_invalid) / len(result.molecule)))
        for cycle_len in sorted(cycle_counts):
            valid = cycle_len == 1 or cycle_len == csm_args['op_order'] or (
                cycle_len == 2 and csm_args['op_type'] == 'SN')
            count = cycle_counts[cycle_len]
            PrintClass.my_print("\tThere %s %d %s %s of length %d" % (
                "is" if count == 1 else "are", count, "invalid" if not valid else "",
                "cycle" if count == 1 else "cycles",
                cycle_len))

    return result


def handle_args(args):
    if not args:
        args = sys.argv[1:]
    if '--approx' not in args:
        args.append("--approx")
    if "--dir" in args:
        raise ValueError("The argument --dir cannot be used in conjunction with test_directions")

    parser = direction_parser()
    parsed_args = parser.parse_args(args)

    direction_choices = [choice_dict[x] for x in parsed_args.direction_choice]
    dirs_file = parsed_args.dirs_file
    k = parsed_args.num_rand
    seed = parsed_args.seed
    if seed:
        seed = int(seed)

    dir_output = parsed_args.dir_output
    PrintClass.set(dir_output)
    args.remove(dir_output)

    args = [x for x in args if x not in parsed_args.direction_choice]
    csm_args = get_split_arguments(args)
    molecule = MoleculeReader.from_file(**csm_args)
    return direction_choices, dirs_file, k, seed, csm_args, molecule


def run(args=[]):
    print("CSM version %s" % __version__)
    direction_choices, dirs_file, k, seed, csm_args, molecule = handle_args(args)

    best_result = CSMState(csm=MAXDOUBLE)
    best_initial_dir = None
    result_csms = []
    result_dirs = []

    for choice in direction_choices:
        dirs = choose_directions(choice, molecule, csm_args, dirs_file, k, seed)

        PrintClass.my_print("Using direction choice", choice, "there are", len(dirs), "initial directions",
                            print_flag=True)
        for dir in dirs:
            result = run_dir(dir, csm_args, molecule)
            if result.csm < best_result.csm:
                best_result = result
                best_initial_dir = dir
            result_csms.append(result.csm)
            result_dirs.append(result.dir)

    PrintClass.my_print("The best csm:", format_CSM(best_result.csm), "\nwas found using the initial dir", best_initial_dir,
                        print_flag=True)
    print_results(best_result, csm_args)
    return result_csms, result_dirs


def run_no_return_dirs(args=[]):
    run(args)


def test_vector_distances():
    def angle_between(v1, v2):
        """
        https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python/13849249#13849249
        """

        def unit_vector(vector):
            """ Returns the unit vector of the vector.  """
            return vector / np.linalg.norm(vector)

        v1_u = unit_vector(v1)
        v2_u = unit_vector(v2)
        return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

    args = [
        # taken from running xyz tests and choosing those with highest volatility
        {
            "args": ['c3',
                     'C:\\Users\\devora.CHELEM\\Sources\\csm\\python\\tests\\xyz_tests\\xyz_mols\\t02_CH4\\molecules\\11.xyz',
                     'nul', '--keep-structure', '--babel-bond'],
            "dir": [-0.94676, -0.068738, 0.314515],
            "csm": 3.4391512605256747
        },
        {
            "args": ['ch',
                     'C:\\Users\\devora.CHELEM\\Sources\\csm\\python\\tests\\xyz_tests\\xyz_mols\\t02_CH4\\molecules\\43.xyz',
                     'nul', '--keep-structure', '--babel-bond'],
            "dir": [-0.239296, 0.480698, -0.843603],
            "csm": 0.04835451722270401
        },
        {
            "args": ['c4',
                     'C:\\Users\\devora.CHELEM\\Sources\\csm\\python\\tests\\xyz_tests\\xyz_mols\\t03_B6H4\\molecules\\43.xyz',
                     'nul', '--keep-structure', '--babel-bond'],
            "dir": [-0.403519, 0.291713, -0.867223],
            "csm": 0.00023937588615741134
        },
        {
            "args": ['c4',
                     'C:\\Users\\devora.CHELEM\\Sources\\csm\\python\\tests\\xyz_tests\\xyz_mols\\t08-cubane\\molecules\\31.xyz',
                     'nul', '--keep-structure', '--babel-bond'],
            "dir": [0.251258, 0.799556, -0.545509],
            "csm": 0.020382748120417737
        },
        {
            'args': ['c3',
                     'C:\\Users\\devora.CHELEM\\Sources\\csm\\python\\tests\\xyz_tests\\xyz_mols\\t08-cubane\\molecules\\31.xyz',
                     'nul', '--keep-structure', '--babel-bond'],
            "dir": [-0.641906, 0.672037, -0.369218],
            "csm": 0.00947757927406423
        }
    ]

    degree = 1
    theta = radians(degree)  # 17
    parts = int(pi * 2 / theta)
    rot_mat_x = np.array([[1, 0, 0],
                          [0, cos(theta), -sin(theta)],
                          [0, sin(theta), cos(theta)]])
    rot_mat_y = np.array([[cos(theta), 0, sin(theta)],
                          [0, 1, 0],
                          [-sin(theta), 0, cos(theta)]])
    rot_mat_z = np.array([[cos(theta), -sin(theta), 0],
                          [sin(theta), cos(theta), 0],
                          [0, 0, 1]])

    for this_arg in args:
        PrintClass.my_print(this_arg)
        csm_args = get_split_arguments(this_arg["args"])
        molecule = MoleculeReader.from_file(**csm_args)
        rotated_dir = this_arg["dir"]
        for i in range(parts):
            PrintClass.my_print("rotated", i, "degrees")
            rotated_dir = rot_mat_z.dot(rotated_dir)
            result = run_dir(rotated_dir, csm_args, molecule)
            if result.csm - this_arg["csm"] > 0.00001:
                PrintClass.my_print("CSM difference:", result.csm - this_arg["csm"])
                PrintClass.my_print("**********", i, "**", angle_between(this_arg["dir"], rotated_dir))
            else:
                PrintClass.my_print("CSM is same")
        PrintClass.my_print("[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]")


if __name__ == '__main__':
    run(args=sys.argv[1:])
    # test_vector_distances()
