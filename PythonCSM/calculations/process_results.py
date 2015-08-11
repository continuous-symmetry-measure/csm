__author__ = 'YAEL'

from calculations.normalizations import normalize_coords,de_normalize_coords


def process_results(results, csm_args):
    """
    Final normalizations and de-normalizations
    :param results: CSM calculations results
    :param csm_args: CSM args
    """
    results.molecule.set_norm_factor(csm_args['molecule'].norm_factor)
    masses = [atom.mass for atom in results.molecule.atoms]
    normalize_coords(results.outAtoms, masses, csm_args['keepCenter'])
    results.molecule.de_normalize()
    results.outAtoms = de_normalize_coords(results.outAtoms, results.molecule.norm_factor)