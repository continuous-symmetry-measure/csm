__author__ = 'YAEL'


class CSMCalculationsData:
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
