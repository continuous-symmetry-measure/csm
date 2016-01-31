__author__ = 'YAEL'


class CSMCalculationsData:
    """ A wrapper class for all the data required during the CSM calculation.
    This class is passed to the old C++ methods.
    To avoid duplication, it's also used by the new Python methods.

    Once the C++ code is gone, this class is going to seem strange. It'll need to be removed.
    """
    # Includes all data passed to/from CPP old_calculations
    def __init__self(self, molecule, op_type, op_order, dir=None, perm=None):
        self.molecule = molecule
        self.outAtoms = []
        self.dir = dir or []
        self.csm = 0
        self.dMin = 0
        self.perm = perm or []
        self.localCSM = []
        self.operationType = op_type
        self.chMinOrder = 2
        self.chMinType = 'CS'
        self.opOrder = op_order
