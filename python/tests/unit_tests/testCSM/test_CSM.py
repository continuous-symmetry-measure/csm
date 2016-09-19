import pytest
import json
from csm.main.csm_run import run
from tests.utils.run_test import close_enough




def test_csm(args, molecule, expected):
    args[1]=molecule
    args[2]=r'D:\UserData\devora\Sources\csm\python\tests\unit_tests\noonecares.txt'
    result=run(args)
    assert close_enough(result.csm, expected.csm)

