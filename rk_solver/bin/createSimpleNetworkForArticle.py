#!/usr/bin/python

from sys import argv
from math import exp

S = float(argv[1])

F1 = 1e-6 * S
F2 = 0.01 * F1
F3 = 0.1 * F1
F4 = 0.01 * F1
F5 = 0.0

W1 = 1e-4
W2 = W1
W3 = 0.1 * W1
W4 = 0.1 * W1
W5 = 0.01 * W1

A1 = 1. / S
A2 = 0.1 * A1
A3 = A1
A4 = 0.01 * A1
A5 = 0.1 * A1

f = file("params.net","w+")
f.write("Begin" + "\n")
f.write("Species\tFlux\tDesorption\tSweepRate\tCutoff" + "\n")
f.write("X1\t" + str(F1) + "\t" + str(W1) + "\t" + str(A1) + "\t" + str(10) + "\n")
f.write("X2\t" + str(F2) + "\t" + str(W2) + "\t" + str(A2) + "\t" + str(10) + "\n")
f.write("X3\t" + str(F3) + "\t" + str(W3) + "\t" + str(A3) + "\t" + str(10) + "\n")
f.write("X4\t" + str(F4) + "\t" + str(W4) + "\t" + str(A4) + "\t" + str(10) + "\n")
f.write("X5\t" + str(F5) + "\t" + str(W5) + "\t" + str(A5) + "\t" + str(10) + "\n")
f.write("Interaction" + "\n")
f.write("X1 X1 => X1-1" + "\n")
f.write("X1 X2 => X1-2" + "\n")
f.write("X1 X3 => X1-3" + "\n")
f.write("X2 X3 => X2-3" + "\n")
f.write("X2 X4 => X3-4" + "\n")
f.write("X3 X5 => X3-5" + "\n")
f.write("End" + "\n")













