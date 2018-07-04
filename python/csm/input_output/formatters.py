import math
import sys

def format_perm_count(num):
    if abs(num) < 1000000:
        return '%d' % num
    return '%.5g' % num

def format_CSM(num):
    return "%.6lf" % (abs(num))

def non_negative_zero(number):
    if math.fabs(number)<0.00001:
        return 0.0000
    else:
        return number

def format_unknown_str(input):
    try:
        num=float(input)
        return "%.4lf" % (num)
    except:
        return str(input)


def csm_log(*strings, file=sys.stderr, **kwargs):
    print(*strings, **kwargs, file=file)