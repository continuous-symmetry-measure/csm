from csm.calculations.approx.approx_calculations import find_best_perm
from csm.calculations.approx.dirs import find_symmetry_directions
from csm.calculations.basic_calculations import process_results
from csm.calculations.constants import MINDOUBLE


def approx_calculation(op_type, op_order, molecule, sn_max=8, use_best_dir=False, get_orthogonal=True, detect_outliers=False,use_chains=False, hungarian=False, print_approx=False, dirs=None, *args, **kwargs):
    if not dirs:
        dirs = find_symmetry_directions(molecule, use_best_dir, get_orthogonal, detect_outliers, op_type)

    if op_type == 'CH':  # Chirality
        #sn_max = op_order
        # First CS
        best_result = find_best_perm('CS', 2, molecule, use_chains, hungarian, print_approx, dirs)
        best_result = best_result._replace(op_type='CS') #unclear why this line isn't redundant
        if best_result.csm > MINDOUBLE:
            # Try the SN's
            for op_order in range(2, sn_max + 1, 2):
                result = find_best_perm('SN', op_order, molecule, use_chains, hungarian, print_approx, dirs)
                if result.csm < best_result.csm:
                    best_result = result._replace(op_type = 'SN', op_order = op_order)
                if best_result.csm < MINDOUBLE:
                    break

    else:
        best_result = find_best_perm(op_type, op_order, molecule, use_chains, hungarian, print_approx, dirs)

    best_result = process_results(best_result)
    #TODO: is calc_local supposed to be supported in approx?
    #if calc_local:
    #    local_csm = compute_local_csm(molecule, best_result.perm, best_result.dir, best_result.op_type,
    #                                  best_result.op_order)
    #    best_result = best_result._replace(local_csm=local_csm)

    return best_result
