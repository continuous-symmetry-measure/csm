__author__ = 'YAEL'


class CSMCalculationsData:
    """ A wrapper class for all the data required during the CSM calculation.
    This class is passed to the old C++ methods.
    To avoid duplication, it's also used by the new Python methods.

    Once the C++ code is gone, this class is going to seem strange. It'll need to be removed.
    """
    # Includes all data passed to/from CPP calculations
    def __init__(self, csm_args=None):
        if csm_args:
            self.molecule = csm_args['molecule']
            self.outAtoms = []
            if 'dir' in csm_args:
                self.dir = csm_args['dir']
            else:
                self.dir = []
            self.csm = 0
            self.dMin = 0
            if 'perm' in csm_args:
                self.perm = csm_args['perm']
            else:
                self.perm = []
            self.localCSM = []
            self.operationType = csm_args['type']
            self.chMinOrder = 2
            self.chMinType = 'CS'
            self.opOrder = csm_args['opOrder']
        else:
            self.molecule = None
            self.outAtoms = []
            self.dir = []
            self.csm = 0
            self.dMin = 0
            self.perm = []
            self.localCSM = []
            self.operationType = ''
            self.chMinOrder = 2
            self.chMinType = 'CS'
            self.opOrder = 0
