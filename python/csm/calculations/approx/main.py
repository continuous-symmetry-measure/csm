"""
This is where the outside wrapper for the approximate calculation is located and called from
"""

import math
from collections import OrderedDict

import datetime
import numpy as np
import sys

from csm.calculations.approx.approximators import HungarianApproximator, GreedyApproximator, ManyChainsApproximator, \
    StructuredApproximator
from csm.calculations.approx.dirs import DirectionChooser
from csm.calculations.data_classes import process_results, Operation, Calculation
from csm.calculations.constants import MINDOUBLE, MAXDOUBLE
from csm.input_output.readers import check_perm_validity
from csm.input_output.formatters import csm_log as print


class ApproxStatistics:
    class DirectionStatistics:
        # per direction, we want to store:
        # 1. every direction passed through
        # 2. every csm passed through, and percent cycle preservation
        # 3. runtime
        def __init__(self, dir, index):
            self.start_dir=dir
            self.index=index
            self.dirs=[]
            self.csms=[]
            self.cycle_stats=[]
            self._stop_reason=""
        def append(self, result):
            self.dirs.append(result.dir)
            self.csms.append(result.csm)
            self.cycle_stats.append(result.num_invalid)

        @property
        def stop_reason(self):
            return self._stop_reason

        @stop_reason.setter
        def stop_reason(self, reason):
            self._stop_reason=reason
        def start_clock(self):
            self.__start_time=datetime.datetime.now()
        def end_clock(self):
            now=datetime.datetime.now()
            time_d = now-self.__start_time
            self.run_time=time_d.total_seconds()
        def __repr__(self):
            return str({
                "dirs":self.dirs,
                "csms":self.csms
            })

        @property
        def end_dir(self):
            return self.dirs[-1]
        @property
        def start_csm(self):
            return self.csms[0]
        @property
        def end_csm(self):
            return self.csms[-1]
        @property
        def num_iterations(self):
            return len(self.dirs)
        def __lt__(self, other):
            return self.end_csm < other.end_csm

    def __init__(self, initial_directions):
        self.directions_dict=OrderedDict()
        self.directions_arr=[]
        for index, dir in enumerate(initial_directions):
            self.directions_dict[tuple(dir)]=ApproxStatistics.DirectionStatistics(dir, index)
            self.directions_arr.append(self.directions_dict[tuple(dir)])
    def __getitem__(self, key):
        return self.directions_dict[tuple(key)]
    def __iter__(self):
        return self.directions_dict.__iter__()
    def __str__(self):
        return str(self.directions_dict)

class ApproxCalculation(Calculation):
    """
    Runs an approximate algorithm to estimate the csm value, using directions to create permutations iteratively
    """
    def __init__(self, operation, molecule, approx_algorithm='hungarian', use_best_dir=False, get_orthogonal=True, detect_outliers=False, dirs=None, keep_structure=False, timeout=100, fibonacci=False, num_dirs=50, selective=False, num_selected=10, *args, **kwargs):
        """
        Initializes the ApproxCalculation
        :param operation: instance of Operation class or named tuple, with fields for name and order, that describes the symmetry
        :param molecule: instance of Molecule class on which the described symmetry calculation will be performed
        :param approx_algorithm: a string, can be 'hungarian', 'greedy', or 'many-chains'. Chooses the approximator alogrithm to use
        :param use_best_dir: use only the best direction 
        :param get_orthogonal: get orthogonal direction vectors from the main directions
        :param detect_outliers: detect outliers and use the imrpoved direction vectors
        :param dirs: a list of directions to use as initial dire (will override all other dir-choice options)
        """
        super().__init__(operation, molecule)
        if approx_algorithm == 'hungarian':
            self.approximator_cls = HungarianApproximator
        if approx_algorithm == 'greedy':
            self.approximator_cls = GreedyApproximator
        if approx_algorithm == 'many-chains':
            self.approximator_cls = ManyChainsApproximator
        if keep_structure:
            self.approximator_cls = StructuredApproximator

        self.timeout=timeout
        dir_chooser=DirectionChooser(molecule, operation.type, operation.order, use_best_dir, get_orthogonal, detect_outliers, dirs, fibonacci, num_dirs)
        self.directions=dir_chooser.dirs
        self.statistics=ApproxStatistics(self.directions)
        self.selective=selective
        self.num_selected=num_selected

    def calculate(self):
        op_type=self.operation.type
        op_order=self.operation.order
        molecule=self.molecule

        # step two: run the appropriate approximator
        if op_type == 'CH':  # Chirality
            # First CS
            approximator = self.approximator_cls('CS', 2, molecule, self.directions, self.statistics, self.log)
            best_result = approximator.approximate()
            if best_result.csm > MINDOUBLE:
                # Try the SN's
                approximator = self.approximator_cls('SN', 2, molecule, self.directions, self.statistics, self.log)
                for op_order in range(2, self.operation.order + 1, 2):
                    approximator._op_order = op_order
                    result = approximator.approximate()
                    if result.csm < best_result.csm:
                        best_result = result
                        best_result.op_type='SN'
                        best_result.op_order=op_order
                    if best_result.csm < MINDOUBLE:
                        break
        else:
            approximator = self.approximator_cls(op_type, op_order, molecule, self.directions, self.statistics, self.log, self.timeout, self.selective, self.num_selected)
            best_result = approximator.approximate()
        # step three: process and return results
        self.approximator=approximator
        self._csm_result= best_result
        return self.result

    def log(self, *args, **kwargs):
        pass


class PrintApprox(ApproxCalculation):
    def log(self, *args, **kwargs):
        print(*args)


