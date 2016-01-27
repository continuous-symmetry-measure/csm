__author__ = 'zmbq'

import math
import logging
logger = logging.getLogger("csm")

def normalize_coords(coords, masses, keep_center):
    """
    Normalize coordinates
    :param coords: atom coordinates
    :param masses: Atomic masses
    :param keep_center:  When false, the center of mass is moved to (0,0,0)
    :return: (List of normalized coordinates, normalization factor)
    """

    x_avg = y_avg = z_avg = 0.0

    if not keep_center:
        mass_sum = 0
        size = len(masses)
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
        #logger.log("Norm: %lf i: %lf temp %lf" % (norm, i, tmp))

    # normalize to 1 and not molecule size
    # norm = sqrt(norm / (double)m->size());

    norm = math.sqrt(norm)
    #logger.log("Second normalization factor is %lf and average is (%lf, %lf, %lf)" % (norm, x_avg, y_avg, z_avg))

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

