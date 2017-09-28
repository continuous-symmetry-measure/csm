import pytest
import json
from csm.main.csm_run import run as csmrun
from tests.utils.run_test import close_enough, YaffaError
import os

class TestClass:
    args=None
    def run_test(self, args, molecule):
        args[1] = molecule
        if self.__class__.args !=list(args):
            self.__class__.args = list(args)
            args[2] = os.devnull
            print(args)

            result = csmrun(args)
            self.__class__.result=result
        return self.__class__.result

class TestCSM(TestClass):
    def test_csm(self, args, molecule, expected, equiv):
        print("CSM")
        result = self.run_test(args, molecule)
        if not close_enough(result.csm, result.formula_csm):
            raise YaffaError
        print(expected.csm)
        assert close_enough(result.csm, expected.csm)

class xTestPermsSame(TestClass):
    def test_perm_same(self, args, molecule, expected, equiv):
        import csv
        args1=list(args)
        args1.append("--output-perms")
        args1.append(r"C:\Users\devora.CHELEM\Sources\olderetired\dumpingground\perms1.csv")
        result = self.run_test(args1, molecule)
        set1=set()

        args2=list(args)
        args2.append("--no-constraint")
        args2.append("--output-perms")
        args2.append(r"C:\Users\devora.CHELEM\Sources\olderetired\dumpingground\perms2.csv")
        result2 = self.run_test(args2, molecule)
        set2=set()

        with open(r"C:\Users\devora.CHELEM\Sources\olderetired\dumpingground\perms1.csv") as csvfile:
            with open(r"C:\Users\devora.CHELEM\Sources\olderetired\dumpingground\perms2.csv") as csvfile2:
                reader = csv.DictReader(csvfile)
                reader2= csv.DictReader(csvfile2)
                for row1, row2 in zip(reader, reader2):
                    set1.add(row1['Permutation'])
                    set2.add(row2['Permutation'])

        assert len(set1)==len(set2)

        set3=set1.union(set2)

        assert len(set3)==len(set2)





class xTestSymm(TestClass):
    def test_symmetric_structure(self, args, molecule, expected, equiv):
        print("SYMM")
        result=self.run_test(args, molecule)
        atom_pos_not_match = []
        for i in range(len(expected.molecule)):
            for index in range(3):
                val1 = expected.symmetric_structure[i][index]
                val2 = result.symmetric_structure[i][index]
                if not close_enough(val1,val2):
                    atom_pos_not_match.append((str(i), index, val1, val2))
        assert(len(atom_pos_not_match)==0)

class xTestPerm(TestClass):
    def test_perm(self, args, molecule, expected, equiv):
        print("PERM")
        result=self.run_test(args, molecule)
        symm=args[0]
        if symm[0].upper()=="S":
            symm = "c" + symm[1:]

        exp_p=expected.perm
        res_p=result.perm
        if expected.perm != result.perm:
            if symm in equiv:
                for pair in equiv[symm]:
                    if expected.perm in pair and result.perm in pair:
                        res_p=exp_p
                        break
        assert(exp_p==res_p)
