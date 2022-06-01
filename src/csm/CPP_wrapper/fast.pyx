# cython: profile=True
# cython: language-level=3
# cython: boundscheck=False, wraparound=False, nonecheck=False

# Include all pyx files, as explained here:  http://stackoverflow.com/a/11804020/871910

include "misc.pxi"
include "cache.pyx"
include "primitives.pyx"
include "fast_calculations.pyx"
include "permuters.pyx"
include "approx.pyx"