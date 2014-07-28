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

def clean_molecules(lines):
    """
    Cleans the molecule output in lines
    Sometimes a molecule looks like:
    ##
    pathname
    At X Y Z
    At X Y Z
    .
    .
    .

    And sometimes the pathname is missing

    This change is due to different OpenBabel versions, and should not affect our tests

    :param lines: Output lines
    :return: Cleaned output lines, without the extra pathname
    """

    output = []
    ignore = None
    for i in range(len(lines)):
        if ignore is not None and ignore==i:
            ignore = None
            continue
        output.append(lines[i])
        if lines[i] in ['INITIAL STRUCTURE COORDINATES', 'RESULTING STRUCTURE COORDINATES']:  # Molecule header
            suspect = lines[i+2]
            if len(suspect.split(' '))==1: # Path, not the first table line
                ignore = i+2

    return output

def compare_files(file1, file2):
    lines1 = read_file(file1)
    lines2 = read_file(file2)

    # Output may not be identical - OpenBabel changed the way molecules are displayed - sometimes the molecule files path
    # is displayed, too. We need to ignore it, as it will make tests look different based on the file-system
    lines1 = clean_molecules(lines1)
    lines2 = clean_molecules(lines2)

    if len(lines1)!=len(lines2):
        return False
    for i in range(len(lines1)):
        if lines1[i]!=lines2[i]:
            return False

    return True

def report_error(expected_file, generated_file):
    print("Files are not identical\n\nEXPECTED:\n")
    for line in read_file(expected_file):
        print (line)
    print ("\n\nGENERATED:")
    for line in read_file(generated_file):
        print (line)

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
