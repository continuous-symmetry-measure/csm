'''
the classes used for running the approximate algorithm
'''

import datetime
import multiprocessing

import numpy as np

from csm.calculations.approx.perm_builders import _OptionalLogger, _HungarianPermBuilder, _GreedyPermBuilder, \
    _ManyChainsPermBuilder, _StructuredPermBuilder
from csm.calculations.approx.statistics import SingleDirectionStatistics, ApproxStatistics
from csm.calculations.basic_calculations import CalculationTimeoutError, check_timeout
from csm.calculations.basic_calculations import now, run_time
from csm.calculations.constants import MAXDOUBLE, CSM_THRESHOLD
from csm.calculations.data_classes import CSMState, CSMResult, BaseCalculation
from csm.calculations.exact_calculations import ExactCalculation

from csm.input_output.formatters import csm_log as print





class SingleDirApproximator(_OptionalLogger):
    def __init__(self, operation, molecule, perm_from_dir_builder, log_func=None, timeout=100,
                 max_iterations=50, chain_perms=None):
        self._log_func = log_func
        self._molecule = molecule
        self._op_type = operation.type
        self._op_order = operation.order
        self._operation = operation
        self.max_iterations = max_iterations
        self.perm_from_dir_builder = perm_from_dir_builder(operation, molecule, log_func, timeout)
        if not chain_perms:
            self._chain_permutations = self.perm_from_dir_builder.get_chain_perms()
        else:
            self._chain_permutations=chain_perms
        self.timeout = timeout
        self.start = datetime.datetime.now()

    def _create_perm_from_dir(self, dir, chainperm):
        return self.perm_from_dir_builder.create_perm_from_dir(dir, chainperm)

    def calculate(self, dir):
        statistics = SingleDirectionStatistics(dir)
        statistics.start_clock()
        best = CSMState(molecule=self._molecule, op_type=self._op_type, op_order=self._op_order, csm=MAXDOUBLE)
        self._log("Calculating for initial direction: ", dir)
        for chainperm in self._chain_permutations:
            if len(self._chain_permutations) > 1:
                self._log("\tCalculating for chain permutation ", chainperm)
            best_for_chain_perm = old_results = CSMState(molecule=self._molecule, op_type=self._op_type,
                                                         op_order=self._op_order,
                                                         csm=MAXDOUBLE, dir=dir)
            i = 0
            while True:
                check_timeout(self.start, self.timeout)
                i += 1

                self._log("\t\titeration", i, ":")

                try:
                    perm = self._create_perm_from_dir(old_results.dir, chainperm)
                    interim_results = ExactCalculation.exact_calculation_for_approx(self._operation,
                                                                                    self._molecule, perm=perm)
                except CalculationTimeoutError as e:
                    self._log("\t\titeration ", i, " Timed out after ", str(e.timeout_delta), " seconds")
                    break

                statistics.append_sub_direction(interim_results)
                # self._log("\t\t\tfound a permutation using dir", old_results.dir, "...")
                if i > 1:
                    self._log("\t\t\tthere are",
                              len(perm) - np.sum(np.array(perm) == np.array(old_results.perm)),
                              "differences between new permutation and previous permutation")
                self._log("\t\t\tthe csm for this permutation is:", str(round(interim_results.csm, 8)))
                self._log("\t\t\tthe new direction from this permutation is", interim_results.dir)
                self._log("\t\t\tthe distance between the new direction and the previous direction is:",
                          str(round(np.linalg.norm(interim_results.dir - old_results.dir), 8)))

                if interim_results.csm < best_for_chain_perm.csm:
                    best_for_chain_perm = interim_results

                # Various stop conditions for the loop, listed as multiple if statements so that the code is clearer
                if i >= self.max_iterations:
                    statistics.stop_reason = "Max iterations"
                    self._log("\t\tStopping after %d iterations" % i)
                    break
                # if i > 1 and math.fabs(old_results.csm - interim_results.csm) / math.fabs(old_results.csm) > 0.01:
                #    self._log("\t\tStopping due to CSM ratio")
                if best_for_chain_perm.csm < CSM_THRESHOLD:
                    statistics.stop_reason = "CSM below threshold"
                    self._log("\t\tStopping because the best CSM is good enough")
                    break
                if abs(np.linalg.norm(interim_results.dir - old_results.dir)) <= 0:
                    statistics.stop_reason = "No change in direction"
                    self._log("\t\tStopping because the direction has not changed")
                    break
                if interim_results.csm >= old_results.csm:  # We found a worse CSM
                    statistics.stop_reason = "No improvement in CSM"
                    self._log("\t\tStopping because CSM did not improve (worse or equal)")
                    break

                old_results = interim_results

            if best_for_chain_perm.csm < best.csm:
                best = best_for_chain_perm
                if best_for_chain_perm.csm < CSM_THRESHOLD:
                    break

        statistics.end_clock()
        return best, statistics


class ApproxCalculation(BaseCalculation, _OptionalLogger):
    def __init__(self, operation, molecule, direction_chooser, approx_algorithm='hungarian',
                 log_func=None, selective=False, num_selected=10, chain_perms=None, *args, **kwargs):

        super().__init__(operation, molecule)

        if log_func==None:
            log_func=self._empty_log
        self._log_func = log_func

        # choose the approximator class
        if approx_algorithm == 'hungarian':
            self.perm_builder = _HungarianPermBuilder
        if approx_algorithm == 'greedy':
            self.perm_builder = _GreedyPermBuilder
        if approx_algorithm == 'many-chains':
            self.perm_builder = _ManyChainsPermBuilder
        if approx_algorithm == 'structured':
            self.perm_builder = _StructuredPermBuilder

        # get the directions
        self._initial_directions = direction_chooser.dirs
        if len(self._initial_directions) > 1:
            self._log("There are", len(self._initial_directions),
                      "initial directions to search for the best permutation")

        self.selective = selective
        self.num_selected = num_selected
        self._single_statistics = ApproxStatistics(self._initial_directions)
        self.statistics={}
        self.chain_perms=chain_perms
        self._max_iterations = 30

    def _empty_log(self, *args):
        return

    def calculate(self, timeout=100, *args, **kwargs):
        self.timeout = timeout
        overall_stats = {}
        best_result = super().calculate(timeout)
        overall_stats["runtime"] = run_time(self.start_time)
        self.result = CSMResult(best_result, self.operation, overall_stats=overall_stats,
                                ongoing_stats={"approx": self.statistics})
        return self.result

    def _calculate(self, operation, timeout):
        if operation.type == 'CI' or (operation.type == 'SN' and operation.order == 2):
            dir = [1.0, 0.0, 0.0]
            if operation.type == 'SN':
                op_msg = 'S2'
            else:
                op_msg = 'CI'
            self._log("Operation %s - using just one direction: %s" % (op_msg, dir))
            best_result = self._calculate_for_directions(operation, [dir], 1)
        else:
            if self.selective:
                self._calculate_for_directions(operation, self._initial_directions, 1)
                best_dirs = []
                sorted_csms = sorted(self._single_statistics.directions_arr)
                for item in sorted_csms[:self.num_selected]:
                    best_dirs.append(item.start_dir)
                    self._log("Running again on the", self.num_selected, "best directions")
                best_result = self._calculate_for_directions(operation, best_dirs, self._max_iterations)

            else:
                best_result = self._calculate_for_directions(operation, self._initial_directions, self._max_iterations)
        self.statistics[operation.op_code]=self._single_statistics.to_dict()
        return best_result

    def _calculate_for_directions(self, operation, dirs, max_iterations):
        best = CSMState(molecule=self.molecule, op_type=operation.type, op_order=operation.order, csm=MAXDOUBLE,
                        num_invalid=MAXDOUBLE)
        single_dir_approximator = SingleDirApproximator(operation, self.molecule,
                                                        self.perm_builder, self._log,
                                                        self.timeout, max_iterations=max_iterations, chain_perms=self.chain_perms)
        for dir in dirs:
            best_result_for_dir, statistics = single_dir_approximator.calculate(dir)
            self._single_statistics[dir] = statistics

            if best_result_for_dir.csm < best.csm:
                best = best_result_for_dir
            elif best_result_for_dir.csm == best.csm and best_result_for_dir.num_invalid<best.num_invalid:
                best = best_result_for_dir

            if best.csm < CSM_THRESHOLD:
                    break
        return best


class ParallelApprox(ApproxCalculation):
    def __init__(self, operation, molecule, direction_chooser,
                 log_func=None, pool_size=0, *args, **kwargs):
        if log_func is not None:
            raise ValueError("Cannot run logging on approx in parallel calculation")
        self.pool_size = pool_size
        if pool_size == 0:
            self.pool_size = multiprocessing.cpu_count() - 1
        super().__init__(operation, molecule, direction_chooser, *args, **kwargs)

    def _calculate(self, operation, timeout):
        if operation.type == 'CI' or (operation.type == 'SN' and operation.order == 2):
            raise ValueError("Please don't use parallel calculation for inversion")
        else:
            if self.selective:
                self.max_iterations = 1
                self._calculate_for_directions(operation, self._initial_directions)
                best_dirs = []
                sorted_csms = sorted(self.statistics.directions_arr)
                for item in sorted_csms[:self.num_selected]:
                    best_dirs.append(item.start_dir)
                self.max_iterations = self._max_iterations
                best_result = self._calculate_for_directions(operation, best_dirs)

            else:
                self.max_iterations = self._max_iterations
                best_result = self._calculate_for_directions(operation, self._initial_directions)

        return best_result

    def _calculate_for_directions(self, operation, dirs):
        pool = multiprocessing.Pool(processes=self.pool_size)
        print("Approximating across {} processes".format(self.pool_size))
        single_dir_approximator = SingleDirApproximator(operation, self.molecule,
                                                        self.perm_builder, self._log,
                                                        self.timeout, max_iterations=self.max_iterations, chain_perms=self.chain_perms)
       
        pool_outputs = pool.map(single_dir_approximator.calculate, dirs)
        pool.close()
        pool.join()
        best_result = CSMState(csm=MAXDOUBLE)
        self.statistics = ApproxStatistics(dirs)
        for (result, statistics) in pool_outputs:
            dir = statistics.start_dir
            self.statistics[dir] = statistics
            if result.csm < best_result.csm:
                best_result = result
        self.result = best_result
        return best_result
