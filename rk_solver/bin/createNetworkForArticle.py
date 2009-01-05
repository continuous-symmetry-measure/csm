#!/usr/bin/python

from sys import argv
from math import exp

# X1 - H
# X2 - O
# X3 - OH
# X4 - H2
# X5 - O2
# X6 - H20
# X7 - CO
# X8 - HCO
# X9 - H2CO
# X10 - H3CO
# X11 - CH3OH
# X12 - CO2

S = float(argv[1])

F1 = 1e-8 * S
F2 = 0.01 * F1
F7 = 0.001 * F1

W1 = 1e-4
W2 = W1
W3 = 0.05 * W1 
W7 = 0.01 * W1
W8 = 0.001 * W1 
W9 = 0.05 * W1
W10 = 0.005 * W1

A1 = 1. / S
A2 = 0.1 * A1 
A3 = 0.01 * A1
A7 = 0.001 * A1
A8 = 0.001 * A1
A9 = 0.0005 * A1
A10 = 0.00005 * A1

f = file("params.net","w+")
f.write("Begin" + "\n")
f.write("Species\tFlux\tDesorption\tSweepRate\tCutoff" + "\n")
f.write("X1\t" + str(F1) + "\t" + str(W1) + "\t" + str(A1) + "\t" + str(10) + "\n")
f.write("X2\t" + str(F2) + "\t" + str(W2) + "\t" + str(A2) + "\t" + str(10) + "\n")
f.write("X3\t" + str(0.0) + "\t" + str(W3) + "\t" + str(A3) + "\t" + str(10) + "\n")
f.write("X7\t" + str(F7) + "\t" + str(W7) + "\t" + str(A7) + "\t" + str(10) + "\n")
f.write("X8\t" + str(0.0) + "\t" + str(W8) + "\t" + str(A8) + "\t" + str(10) + "\n")
f.write("X9\t" + str(0.0) + "\t" + str(A9) + "\t" + str(A9) + "\t" + str(10) + "\n")
f.write("X10\t" + str(0.0) + "\t" + str(A10) + "\t" + str(A10) + "\t" + str(10) + "\n")
f.write("Interaction" + "\n")
f.write("X1 X1 => X4" + "\n")
f.write("X1 X2 => X3" + "\n")
f.write("X2 X2 => X5" + "\n")
f.write("X1 X3 => X6" + "\n")
f.write("X1 X7 => X8" + "\n")
f.write("X1 X8 => X9" + "\n")
f.write("X1 X9 => X10" + "\n")
f.write("X1 X10 => X11" + "\n")
f.write("X2 X8  => X12 X1" + "\n")
f.write("X2 X7  => X12" + "\n")
f.write("End" + "\n")













