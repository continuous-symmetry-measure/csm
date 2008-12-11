#!/usr/bin/python

from os import system

st = 10 ** 0.5
grain_sizes = [100, 100 * st, 1000, 1000 * st, 10000, 10000 * st, 100000, 100000 * st]
filenames = []
f = open("sizes", "w")

i = 1
for gs in grain_sizes:
	f.write(str(i) + " " + str(gs) + "\n")
	system("./createMethanolNetwork.py " + str(gs) +"; mv params.net params" + str(i) + ".net" )
	filenames.append("params" + str(i) + ".net")
	i = i+1

f.close()

# run pairs
i = 1
for ff in filenames:
	system("./pairs " + ff + " > pair" + str(i) + ".res; mv pairs.out pairs" + str(i) + ".out")
	i = i+1

# run master
i = 1
for ff in filenames:
	system("./master " + ff + " > master" + str(i) + ".res; mv master.out master" + str(i) + ".out")
	i = i+1

