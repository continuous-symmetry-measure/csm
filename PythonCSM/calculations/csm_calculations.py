__author__ = 'YAEL'

from CPP_wrapper import csm


def perform_operation(csm_args, data):
    if 'dir' in csm_args:
        result = csm.FindBestPermUsingDir(data)
    else:
        if csm_args['findPerm']:
            result = csm.FindBestPerm(data)
        else:
            result = csm.CsmOperation(data)
    return result
