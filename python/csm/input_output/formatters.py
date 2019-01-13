import math
import sys


def format_perm_count(num):
    if abs(num) < 1000000:
        return '%d' % num
    return '%.5g' % num


def format_CSM(num):
    return "%.4lf" % (abs(num))


def non_negative_zero(number):
    if math.fabs(number) < 0.00001:
        return 0.0000
    else:
        return number


def format_unknown_str(input):
    try:
        num = float(input)
        if int(input) == num:
            return int(input)
        return "%.4lf" % (num)
    except:
        return str(input)


global csm_out_pipe
csm_out_pipe = sys.stdout


class output_strings:
    vec = []
    index = 0
    silent = False

    @staticmethod
    def append(s):
        output_strings.vec.append(s)

    @staticmethod
    def fetch():
        if output_strings.index >= len(output_strings.vec):
            return None
        ret = output_strings.vec[output_strings.index]
        output_strings.index += 1
        return ret

    @staticmethod
    def silent_print(*strings, file=csm_out_pipe, **kwargs):
        output_strings.append(" ".join(str(s) for s in strings))
        if not output_strings.silent:
            print(*strings, **kwargs, file=file)


def csm_log(*strings, file=csm_out_pipe, **kwargs):
    output_strings.append(" ".join(str(s) for s in strings))
    print(*strings, **kwargs, file=file)


def silent_print(*strings, file=csm_out_pipe, **kwargs):
    output_strings.silent_print(*strings, **kwargs, file=file)
