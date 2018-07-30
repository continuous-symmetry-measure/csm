# Constants used in the various calculations
import datetime

MINDOUBLE = 1e-8
MAXDOUBLE = 100000000.0
ZERO_IM_PART_MAX = 1e-3
CSM_THRESHOLD=0.0001

class CalculationTimeoutError(TimeoutError):
    def __init__(self, timeout_delta, *args, **kwargs):
        super().__init__("Calculation timed out after "+str(timeout_delta)+" seconds", *args, **kwargs)
        self.timeout_delta=timeout_delta
