import pprint
from tempfile import NamedTemporaryFile, mkstemp

__author__ = 'zmbq'

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

def get_test_args(test_dir, output_path):
    """
    Returns the arg list as should be passed to the subprocess functions
    :param test_dir: Folder with the test files
    :param output_path: The path for the CSM output file
    :return: An argument list, including the CSM executable
    """
    input_args = get_first_line(os.path.join(test_dir, 'args.txt'))
    # args.txt file looks like:
    # mode input_file [-option -option -option...]
    #
    # The actual arguments to CSM should be
    # mode input_file output_file [-option -option -option]
    input_args_list = input_args.split(' ')

    # Parse the first input arguments
    mode = input_args_list[0]
    input_file = input_args_list[1]
    input_path = os.path.join(test_dir, input_file)

    # Make sure input_args[2:] all begin with -
    for i in range(2,len(input_args_list)):
        if input_args_list[i][0]!='-':
            raise ValueError("Unexpected argument, options should be preceeded by a -")

    # Now return the modified arguments
    return [config.CSM_PATH, mode, input_path, output_path] + input_args_list[2:]

def run_csm(test_dir, output_path):
    args = get_test_args(test_dir, output_path)
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
    result['symmetry'] = lines[0]
    result['scaling'] = lines[1]
    result['direction'] = ()
    for i in range(len(lines)):
        line = lines[i]
        if line=='DIRECTIONAL COSINES:':
            result['direction'] = parse_vector(lines[i+1])
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
    if not same_vector(output1['direction'], output2['direction']):
        opposite = tuple([-x for x in output1['direction']])
        if not same_vector(opposite, output2['direction']):
            return False

    return True

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
