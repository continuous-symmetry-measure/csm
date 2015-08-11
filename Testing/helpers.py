import pprint
import sys
from tempfile import NamedTemporaryFile, mkstemp

__author__ = 'zmbq'
import csm

"""
Helper functions for running CSM tests
"""

import config
import os.path
import os
import subprocess

def read_file(path, include_empty=False):
    """
    Returns a list of all lines in the file.
    :param path: File path
    :param include_empty: Whether to return empty lines
    :return: A list of all lines in the file. Lines are stripped for leading and trailing whitespace
    """
    with open(path, "r") as file:
        lines = file.readlines()
    stripped = [line.strip() for line in lines]
    if include_empty:
        return stripped
    return [line for line in stripped if line]

def get_first_line(path):
    """
    Returns the first non empty line of a file
    :param path: File path
    :return: The first non empty line
    """
    return read_file(path)[0]

def get_python_test_args(test_dir, output_path):
    input_args = get_first_line(os.path.join(test_dir, 'args.txt'))
    input_args_list = input_args.split(' ')
    mode = input_args_list[0]
    input_file = input_args_list[1]
    input_path = os.path.join(test_dir, input_file)

    # The Python CSM requires options to start with '--', while the test cases start with '-'
    for i in range(2,len(input_args_list)):
        if input_args_list[i][0] == '-':
            if input_args_list[i][1]!='-':
                input_args_list[i] = '-' + input_args_list[i]  # Prepend another -
        else:
            # This can be a filename, we need to add the path.
            if input_args_list[i-1] in ('--useperm', '--usedir', '--log'):
                input_args_list[i] = os.path.join(test_dir, input_args_list[i])

    # TODO: Add the path for dirfile, permfile and log
    return [mode, input_path, output_path] + input_args_list[2:]

def run_csm(test_dir, output_path):
    args = get_python_test_args(test_dir, output_path)
    print('Executing %s...' % args)
    with subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as ps:
        output = ps.stdout.readlines()
        error = ps.stderr.readlines()
    return (output, error)

def parse_vector(line):
    """
    Parses a 3D vector, returing (x,y,z)
    :param line: A line with a vector of 3 floating point numbers
    :return: (x,y,z)
    """
    floats = [float(s) for s in line.split()]
    return tuple(floats)

def parse_output(lines):
    """
    Parses the output, returning a dictionary with

    symmetry: The symmetry measure line (string)
    scaling: The scaling factor line (string)
    direction: The directional cosine vector (x,y,z)
    :param lines: File lines
    :return: A dictionary with the parsed information
    """
    result = {}
    result['symmetry'] = float(lines[0].split()[-1])
    result['scaling'] = float(lines[1].split()[-1])
    result['direction'] = ()
    result['permutation'] = ()
    for i in range(len(lines)):
        line = lines[i]
        if line=='DIRECTIONAL COSINES:':
            result['direction'] = parse_vector(lines[i+1])
        if line=='PERMUTATION:':
            result['permutation'] = parse_vector(lines[i+1])
    return result


def compare_files(file1, file2):
    def same_vector(vec1, vec2):
        if len(vec1)!=len(vec2):
            return False
        for i in range(len(vec1)):
            num1 = vec1[i]
            num2 = vec2[i]
            if abs(num1-num2) > 1e-6:
                return False
        return True

    output1 = parse_output(read_file(file1))
    output2 = parse_output(read_file(file2))

    if output1['symmetry']!=output2['symmetry']:
        return False
    if output1['scaling']!=output2['scaling']:
        return False
    if not same_vector(output1['permutation'], output2['permutation']):
        return False
    #if not same_vector(output1['direction'], output2['direction']):
    #    opposite = tuple([-x for x in output1['direction']])
    #    if not same_vector(opposite, output2['direction']):
    #        return False

    return True

def compare_results(expected_filename, results):
    expected = parse_output(read_file(expected_filename))

    ok = True
    # Check symmetry
    if abs(expected['symmetry'] - results.csm) > 1e-4:
        print("Expected CSM of %f, got %f" % (expected['symmetry'], results.csm), file=sys.stderr)
        ok = False

    # Check scale
    if abs(expected['scaling'] - results.dMin) > 1e-4:
        print("Expected scaling of %f, got %f" % (expected['scaling'], results.dMin), file=sys.stderr)
        ok = False

    # Check permutation
    expected_perm = [int(p)-1 for p in expected['permutation']]
    if expected_perm!=results.perm:
        print("Expected permutation %s, got %s" % (expected_perm, results.perm), file=sys.stderr)
        ok = False

    # Check direction, first compare straight on
    diff = 0.0
    diff_neg = 0.0
    for i in range(3):
        diff += (expected['direction'][i] - results.dir[i]) ** 2
        diff_neg += (expected['direction'][i] + results.dir[i]) ** 2
    if diff > 1e-5 and diff_neg > 1e-5: # Try inversed
        print("Expected direction %s, got %s" % (expected['direction'], results.dir), file=sys.stderr)
        ok = False

    return ok


def report_error(expected_file, generated_file):
    print("Output was not identical")
    expected = parse_output(read_file(expected_file))
    generated = parse_output(read_file(generated_file))
    pp = pprint.PrettyPrinter(indent=4)
    print("Expected:")
    pp.pprint(expected)
    print("\nGenerated:")
    pp.pprint(generated)

def run_test(test_dir):
    try:
        (tmp_fd, tmp_name) = mkstemp(text=True)
        os.close(tmp_fd)
        (output, error) = run_csm(test_dir, tmp_name)
        ok = compare_files(os.path.join(test_dir, 'result.txt'), tmp_name)
        if not ok:
            report_error(os.path.join(test_dir, 'result.txt'), tmp_name)
            return False
        else:
            print('OK')
            return True
    finally:
        os.remove(tmp_name)

def run_csm_python(test_dir, output_path):
    args = get_python_test_args(test_dir, output_path)
    return csm.run_csm(args)

def run_test_python(test_dir):
    try:
        (tmp_fd, tmp_name) = mkstemp(text=True)
        os.close(tmp_fd)
        results = run_csm_python(test_dir, tmp_name)
        return compare_results(os.path.join(test_dir, 'result.txt'), results)
    finally:
        # os.remove(tmp_name)
        pass
