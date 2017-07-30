from os import path
import pytest
from csm.main.normcsm import normrun
from tests.utils.run_test import close_enough
from conftest import test_folder, output_file, my_tolerance


#os.chdir(test_folder)

@pytest.mark.parametrize("normalization, expected_csm, expected_norm",
                         [
                             (0, 0.0003636, 577.5266368),
                             (1, 0.0110095, 19.0744003),
                             (2, 0.0110084, 19.0763578),
                             (3, 0.0110084, 19.0763578),
                             (4, 0.0110084, 19.0763513),
                             (5, 1.3019489, 12.0000000),
                             (6, 0.1878611, 83.1645538)
                         ])
def test_1iieg3(normalization, expected_csm, expected_norm):
    normalization=str(normalization)
    args=[normalization, 'c3', path.join(test_folder, '1iie-g3.pdb'), output_file, '--approx', '--use-chains']
    result=normrun(args)
    result_norm, result_csm =  result[normalization]
    assert close_enough(result_norm, expected_norm, tolerance=my_tolerance)
    assert close_enough(result_csm, expected_csm, tolerance=my_tolerance)

@pytest.mark.parametrize("normalization, expected_csm, expected_norm",
                         [
                             (0, 0.0096673, 15239.807),
                             (1, 2.8447955, 51.7887132),
                             (2, 2.7667493, 53.2496031),
                             (3, 2.7667493, 53.2496031),
                             (4, 2.7667227, 53.2501147),
                             (5, 28.5790763, 18.0),
                             (6, 0.9824278, 523.6246111)
                         ])
def test_2rlas3(normalization, expected_csm, expected_norm):
    normalization=str(normalization)
    args=[normalization, 'c3', path.join(test_folder, '2rla-s3.pdb'), output_file, '--approx', '--use-chains']
    result=normrun(args)
    result_norm, result_csm =  result[normalization]
    assert close_enough(result_norm, expected_norm, tolerance=my_tolerance)
    assert close_enough(result_csm, expected_csm, tolerance=my_tolerance)

@pytest.mark.parametrize("normalization, expected_csm, expected_norm",
                         [
                             (0, 0.0062539,413.710),
                             (1, 0.0169407, 152.7268056),
                             (2, 0.0169380, 152.7508759),
                             (3, 0.0169380, 152.7508759),
                             (4, 0.0169380, 152.7513621),
                             (5, 3.0487503, 27.0),
                             (6, 0.7931148, 103.7885845)
                         ])
def test_2m7wq3(normalization, expected_csm, expected_norm):
    normalization=str(normalization)
    args=[normalization, 'c3', path.join(test_folder, '2m7w-q3.pdb'), output_file, '--approx', '--use-chains']
    result=normrun(args)
    result_norm, result_csm =  result[normalization]
    assert close_enough(result_norm, expected_norm, tolerance=my_tolerance)
    assert close_enough(result_csm, expected_csm, tolerance=my_tolerance)