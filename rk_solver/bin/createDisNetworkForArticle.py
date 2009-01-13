#!/usr/bin/python

from sys import argv
from math import exp

#X1 = H
#X2 = C
#X3 = CH
#X4 = C2H
#X5 = CH2

S = float(argv[1])

F1 = 1e-7 * S
F2 = 0.5 * F1
F3 = 0.1 * F1

W1 = 1e-5
W2 = W1
W3 = 0.1 * W1

A1 = 0.5 / S
A2 = 0.1 * A1
A3 = A1

D = 5e-4

f = file("params.net","w+")
f.write("Begin" + "\n")
f.write("Species\tFlux\tDesorption\tSweepRate\tCutoff" + "\n")
f.write("X1\t" + str(F1) + "\t" + str(W1) + "\t" + str(A1) + "\t" + str(10) + "\n")
f.write("X2\t" + str(F2) + "\t" + str(W2) + "\t" + str(A2) + "\t" + str(10) + "\n")
f.write("X3\t" + str(F3) + "\t" + str(W3) + "\t" + str(A3) + "\t" + str(10) + "\n")
f.write("Interaction" + "\n")
f.write("X1 X2 => X3" + "\n")
f.write("X1 X3 => X5" + "\n")
f.write("X2 X3 => X4" + "\n")
f.write("Dissociation" + "\n")
f.write("X3 => X1 X2 | " + str(D) + "\n")
f.write("End" + "\n")













