#!/usr/bin/python

from os import system

st = 10 ** 0.5
grain_sizes = [10, 10 * st, 100, 100 * st, 1000, 1000 * st, 10000, 10000 * st, 100000]
filenames = []
f = open("sizes", "w")

i = 1
for gs in grain_sizes:
	f.write(str(i) + " " + str(gs) + "\n")
	system("bash ./createSimpleNetworkForArticle.py " + str(gs) +"; mv params.net params" + str(i) + ".net" )
	filenames.append("params" + str(i) + ".net")
	i = i+1

f.close()

# run rate
i = 1
for ff in filenames:
	system("ulimit -t 200; ./rate " + ff + " > rate" + str(i) + ".res; mv rate.out rate" + str(i) + ".out")
	i = i+1

# run pairs
i = 1
for ff in filenames:
	system("ulimit -t 2000; ./pairs " + ff + " > pair" + str(i) + ".res; mv pairs.out pairs" + str(i) + ".out")
	i = i+1

# run master
i = 1
for ff in filenames:
	system("ulimit -t 100000; ./master " + ff + " > master" + str(i) + ".res; mv master.out master" + str(i) + ".out")
	i = i+1

