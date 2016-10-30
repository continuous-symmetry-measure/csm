from csm.calculations.constants import MINDOUBLE

__author__ = 'zmbq'

import math
import logging
logger = logging.getLogger("csm")

def normalize_coords(coords, masses):
    """
    Normalize coordinates
    :param coords: atom coordinates
    :param masses: Atomic masses
    :return: (List of normalized coordinates, normalization factor)
    """

    x_avg = y_avg = z_avg = 0.0

    size = len(masses)

    mass_sum = 0
    for i in range(size):
        x_avg += coords[i][0] * masses[i]
        y_avg += coords[i][1] * masses[i]
        z_avg += coords[i][2] * masses[i]
        mass_sum += masses[i]
    x_avg /= mass_sum
    y_avg /= mass_sum
    z_avg /= mass_sum

    norm = 0.0
    for i in range(size):
        tmp = (coords[i][0] - x_avg) ** 2 + (coords[i][1] - y_avg) ** 2 + (coords[i][2] - z_avg) ** 2
        norm += tmp
        logger.debug("Norm: %lf i: %lf temp %lf" % (norm, i, tmp))

    # normalize to 1 and not molecule size
    # norm = sqrt(norm / (double)m->size());

    norm = math.sqrt(norm)
    logger.debug("Second normalization factor is %lf and average is (%lf, %lf, %lf)" % (norm, x_avg, y_avg, z_avg))

    if norm<=MINDOUBLE: #in the original code, this check was against MINDOUBLE.
        raise(ValueError("Normalization factor equals zero"))
        #norm=default_value


    for i in range(size):
        coords[i] = ((coords[i][0] - x_avg) / norm, (coords[i][1] - y_avg) / norm, (coords[i][2] - z_avg) / norm)


    return coords, norm


def de_normalize_coords(coords, norm_factor):
    size = len(coords)
    for i in range(size):
        coords[i] = (coords[i][0] * norm_factor,
                     coords[i][1] * norm_factor,
                     coords[i][2] * norm_factor)
    return coords

