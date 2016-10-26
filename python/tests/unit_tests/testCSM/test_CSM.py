import pytest
import json
from csm.main.csm_run import run as csmrun
from tests.utils.run_test import close_enough


class TestClass:
    args=None
    def run_test(self, args, molecule):
        print(args)
        args[1] = molecule
        if self.__class__.args !=list(args):
            self.__class__.args = list(args)
            args[2] = r'D:\UserData\devora\Sources\csm\python\tests\unit_tests\noonecares.txt'
            result = csmrun(args)
            self.__class__.result=result
        return self.__class__.result

    def test_csm(self, args, molecule, expected, equiv):
        print("CSM")
        result = self.run_test(args, molecule)
        print("res",result.csm, "exp", expected.csm)
        assert close_enough(result.csm, expected.csm)

    def test_symmetric_structure(self, args, molecule, expected, equiv):
        print("SYMM")
        result=self.run_test(args, molecule)
        atom_pos_not_match = []
        for i in range(len(expected.molecule)):
            for index in range(3):
                val1 = expected.symmetric_structure[i][index]
                val2 = result.symmetric_structure[i][index]
                print("exp",val1,"res", val2)
                if not close_enough(val1,val2):
                    atom_pos_not_match.append((str(i), index, val1, val2))
        assert(len(atom_pos_not_match)==0)

    def test_perm(self, args, molecule, expected, equiv):
        print("PERM")
        result=self.run_test(args, molecule)
        symm=args[0]
        if symm[0].upper()=="S":
            symm = "c" + symm[1:]

        exp_p=expected.perm
        res_p=result.perm
        print("res", res_p, "exp", exp_p)
        if expected.perm != result.perm:
            if symm in equiv:
                for pair in equiv[symm]:
                    if expected.perm in pair and result.perm in pair:
                        res_p=exp_p
                        break
        assert(exp_p==res_p)
