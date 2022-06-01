import datetime
from collections import OrderedDict

from csm.calculations.constants import MAXDOUBLE
from csm.calculations.data_classes import CSMState, get_chain_perm_string


class SingleDirectionStatistics:
    # per direction, we want to store:
    # 1. every direction passed through
    # 2. every csm passed through, and percent cycle preservation
    # 3. runtime
    def __init__(self, dir):
        self.start_dir = dir
        self.results=[]

        self.dirs = [] #kept for tests to prove algorithms different
        self._stop_reason = ""
        self.least_invalid = CSMState(csm=MAXDOUBLE, num_invalid=MAXDOUBLE)

    def append_sub_direction(self, result):
        self.results.append(result)
        self.dirs.append(result.dir)
        if result.num_invalid < self.least_invalid.num_invalid or \
                (result.num_invalid == self.least_invalid.num_invalid
                 and result.csm < self.least_invalid.csm):
            self.least_invalid = result

    @property
    def stop_reason(self):
        return self._stop_reason

    @stop_reason.setter
    def stop_reason(self, reason):
        self._stop_reason = reason

    def start_clock(self):
        self.__start_time = datetime.datetime.now()

    def end_clock(self):
        now = datetime.datetime.now()
        time_d = now - self.__start_time
        self.run_time = time_d.total_seconds()

    def __repr__(self):
        return "stats for dir"+str(self.start_dir)

    @property
    def end_dir(self):
        return self.results[-1].dir

    @property
    def start_csm(self):
        return self.results[0].csm

    @property
    def end_csm(self):
        return self.results[-1].csm

    @property
    def num_iterations(self):
        return len(self.results)

    @property
    def chain_perm(self):
        arr, st= get_chain_perm_string(self.results[-1].molecule, self.results[-1].perm)
        return st

    @property
    def validity_dict(self):
        valid_struc= (1 - (self.results[-1].num_invalid / len(self.least_invalid.molecule))) * 100
        best_valid=(1 - (self.least_invalid.num_invalid / len(self.least_invalid.molecule))) * 100

        if best_valid > valid_struc:
            return {
            "res valid": valid_struc,
            "dir": self.least_invalid.dir,
            "csm": self.least_invalid.csm,
            "per": best_valid
            }
        else:
            return {
                "res valid": valid_struc,
                "dir": "n/a",
                "csm": "n/a",
                "per": "n/a"
            }


    def __lt__(self, other):
        try:  # if other doesn't have csm, we are less than them
            that_one = other.end_csm
        except:
            return True

        try:  # if we don't have csm, other is less than us
            this_one = self.end_csm
        except:
            return False

        return self.end_csm < other.end_csm

    def to_dict(self):
        try:
            return_dict = {
                "start dir": list(self.start_dir),
                "start csm": self.start_csm,
                "stop reason": self.stop_reason,
                "end dir": list(self.end_dir),
                "end csm": self.end_csm,
                "num iterations": self.num_iterations,
                "run time": self.run_time,
                "chain perm":self.chain_perm,
                "validity":self.validity_dict,
                "dirs":self.dirs
            }
            return return_dict
        except:
            return {
                "start dir": list(self.start_dir),
                "stop reason": "was never reached"
            }


class DirectionStatisticsContainer:
    def __init__(self, initial_directions):
        self.directions_dict = OrderedDict()
        self.directions_arr = []
        for index, dir in enumerate(initial_directions):
            self.directions_dict[tuple(dir)] = SingleDirectionStatistics(dir)
            self.directions_arr.append(self.directions_dict[tuple(dir)])

    def __getitem__(self, key):
        return self.directions_dict[tuple(key)]

    def __setitem__(self, key, value):
        self.directions_dict[tuple(key)] = value

    def __iter__(self):
        return self.directions_dict.__iter__()

    def __str__(self):
        return str(self.directions_dict)

    def to_dict(self):
        return [{"dir": dir, "stats": self.directions_dict[dir].to_dict()} for dir in self.directions_dict]



class ApproxStatistics(DirectionStatisticsContainer):
    pass
