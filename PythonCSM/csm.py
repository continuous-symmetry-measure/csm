"""
Performs some tests on the CSM C++ wrapper
"""
__author__ = 'zmbq'

from CPP_wrapper import csm
from arguments import process_arguments, create_parser

if __name__=='__main__':
    parser = create_parser()
    result = parser.parse_args() # Parse sys.args
    args = process_arguments(result)

    csm.RunCSM(args)
