#!/usr/bin/python

from sys import argv
from math import exp

d = float(argv[1]) / 10000
s = 5e13 # site density per cm^2
S = (3.1415926535897931 * d * d * s)

T = 15

F1 = 2.75e-9 * S
F2 = 0.01 * F1
F7 = 0.001 * F1

NU = 1e12 
k = 8.617343e-5

W1 = NU * (exp(-52e-3 / k / T))
W2 = NU * (exp(-54e-3 / k / T))
W3 = NU * (exp(-54e-3 / k / T))

A1 = NU * (exp(-44e-3 / k / T) / S)
A2 = NU * (exp(-47e-3 / k / T) / S)
A3 = NU * (exp(-47e-3 / k / T) / S)

f = file("params.net","w+")
f.write("Begin" + "\n")
f.write("Species\tFlux\tDiffusion\tSweepRate\tCutoff" + "\n")
f.write("H\t" + str(F1) + "\t" + str(W1) + "\t" + str(A1) + "\t" + str(10) + "\n")
f.write("O\t" + str(F2) + "\t" + str(W2) + "\t" + str(A2) + "\t" + str(10) + "\n")
f.write("OH\t" + str(0.0) + "\t" + str(W3) + "\t" + str(A3) + "\t" + str(10) + "\n")
f.write("Interaction" + "\n")
f.write("H H => H2" + "\n")
f.write("H O => OH" + "\n")
f.write("O O => O2" + "\n")
f.write("H OH => H2O" + "\n")
f.write("End" + "\n")













