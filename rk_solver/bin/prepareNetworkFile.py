#!/usr/bin/python

from sys import argv
from math import exp

S = float(argv[1])
T = 10

F1 = 2.75e-9 * S
F2 = 0.01 * F1
F7 = 0.001 * F1

NU = 1e12 
k = 8.617343e-5
W1 = NU * (exp(-32e-3 / k / T))
W2 = NU * (exp(-32e-3 / k / T))
W3 = NU * (exp(-35e-3 / k / T))
W7 = NU * (exp(-36e-3 / k / T))
W8 = NU * (exp(-38e-3 / k / T))
W9 = NU * (exp(-39e-3 / k / T))
W10 = NU * (exp(-41e-3 / k / T))

A1 = NU * (exp(-22e-3 / k / T) / S)
A2 = NU * (exp(-25e-3 / k / T) / S)
A3 = NU * (exp(-28e-3 / k / T) / S)
A7 = NU * (exp(-30e-3 / k / T) / S)
A8 = NU * (exp(-30e-3 / k / T) / S)
A9 = NU * (exp(-33e-3 / k / T) / S)
A10 = NU * (exp(-35e-3 / k / T) / S)

f = file("params.net","w+")
f.write("Begin" + "\n")
f.write("Species\tFlux\tDiffusion\tSweepRate\tCutoff" + "\n")
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













