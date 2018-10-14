from csm.calculations.constants import MINDOUBLE

__author__ = 'zmbq'

import math
import logging
logger = logging.getLogger("csm")


def calculate_norm_factor(coords, center_of_mass):
    """
    :param coords: list of atom coordinates
    :param center_of_mass: center of mass of molecule
    """
    size=len(coords)
    norm = 0.0
    for i in range(size):
        tmp = (coords[i][0] - center_of_mass[0]) ** 2 + (coords[i][1] - center_of_mass[1]) ** 2 + (coords[i][2] - center_of_mass[2]) ** 2
        norm += tmp
        #logger.debug("Norm: %lf i: %lf temp %lf" % (norm, i, tmp))

    # normalize to 1 and not molecule size
    # norm = sqrt(norm / (double)m->size());

    norm = math.sqrt(norm)
    #logger.debug("Second normalization factor is %lf and average is (%lf, %lf, %lf)" % (norm, center_of_mass[0], center_of_mass[1], center_of_mass[2]))

    if norm<=MINDOUBLE: #in the original code, this check was against MINDOUBLE.
        raise(ValueError("Normalization factor equals zero"))
        #norm=default_value

    return norm

def normalize_coords(coords, center_of_mass, norm_factor):
    """
    :param coords: list of atom coordinates
    :param center_of_mass: center of mass of molecule
    :param norm_factor: nromalization factor (can be calculated with calculate_norm_factor)
    :return: list of atom coordinates, normalized and moved to origin
    """

    for i in range(len(coords)):
        coords[i] = ((coords[i][0] - center_of_mass[0]) / norm_factor,
                     (coords[i][1] - center_of_mass[1]) / norm_factor,
                     (coords[i][2] - center_of_mass[2]) / norm_factor)


    return coords


def de_normalize_coords(coords, norm_factor):
    size = len(coords)
    for i in range(size):
        coords[i] = (coords[i][0] * norm_factor,
                     coords[i][1] * norm_factor,
                     coords[i][2] * norm_factor)
    return coords

