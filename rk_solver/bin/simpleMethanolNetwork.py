#!/usr/bin/python

from sys import argv
from math import exp

base_flux = float(argv[1])
baseA = float(argv[2]);


F1 = base_flux
F2 = 0.01 * F1
F7 = 0.001 * F1

W1 = 1
W2 = 1
W3 = 1
W7 = 1
W8 = 1
W9 = 1
W10 = 1

A1 = baseA
A2 = 0.1 * baseA
A3 = 0.1 * baseA
A7 = 0.1 * baseA
A8 = 0.1 * baseA
A9 = 0.1 * baseA
A10 = 0.1 * baseA

f = file("params.net","w+")
f.write("Begin" + "\n")
f.write("Species\tFlux\tDesorption\tSweepRate\tCutoff" + "\n")
f.write("H\t" + str(F1) + "\t" + str(W1) + "\t" + str(A1) + "\t" + str(10) + "\n")
f.write("O\t" + str(F2) + "\t" + str(W2) + "\t" + str(A2) + "\t" + str(10) + "\n")
f.write("OH\t" + str(0.0) + "\t" + str(W3) + "\t" + str(A3) + "\t" + str(10) + "\n")
f.write("CO\t" + str(F7) + "\t" + str(W7) + "\t" + str(A7) + "\t" + str(10) + "\n")
f.write("HCO\t" + str(0.0) + "\t" + str(W8) + "\t" + str(A8) + "\t" + str(10) + "\n")
f.write("H2CO\t" + str(0.0) + "\t" + str(A9) + "\t" + str(A9) + "\t" + str(10) + "\n")
f.write("H3CO\t" + str(0.0) + "\t" + str(A10) + "\t" + str(A10) + "\t" + str(10) + "\n")
f.write("Interaction" + "\n")
f.write("H H => H2" + "\n")
f.write("H O => OH" + "\n")
f.write("O O => O2" + "\n")
f.write("H OH => H2O" + "\n")
f.write("H CO => HCO" + "\n")
f.write("H HCO => H2CO" + "\n")
f.write("H H2CO => H3CO" + "\n")
f.write("H H3CO => CH3OH" + "\n")
f.write("End" + "\n")













