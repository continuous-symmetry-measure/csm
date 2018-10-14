# Constants used in the various calculations
from datetime import datetime

MINDOUBLE = 1e-8
MAXDOUBLE = 100000000.0
ZERO_IM_PART_MAX = 1e-3
CSM_THRESHOLD=0.0001

global global_start_time
global_start_time= datetime.now()

global global_time_out
global_time_out=50000


def set_global_timeout(timeout):
    global global_time_out
    global_time_out=timeout