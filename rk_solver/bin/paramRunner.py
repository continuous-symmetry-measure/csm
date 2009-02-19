#!/usr/bin/python

from os import system
from sys import argv

st = 10 ** 0.5
grain_sizes = [10, 10 * st, 100, 100 * st, 1000, 1000 * st, 10000, 10000 * st, 100000]
indices = [1,2,3,4,5,6,7]


filenames = []
f = open("sizes", "w")

i = 1
for gs in grain_sizes:
	f.write(str(i) + " " + str(gs) + "\n")
	system("./" + argv[1] + " " + str(gs) +"; mv params.net params" + str(i) + ".net" )
	filenames.append("params" + str(i) + ".net")
	i = i+1

f.close()

# run rate
system("rm -f rate.dat")
for i in indices:	
	ff = filenames[i - 1]
	system("ulimit -t 200; ./rate " + ff + " > rate" + str(i) + ".res; mv rate.out rate" + str(i) + ".out")
	system("echo " + str(grain_sizes[i]) + " `tail -1 rate" + str(i) + ".out | cut -f3-` >> rate.dat") 

# run master
system("rm -f master.dat")
master_cmd = "master_iter"
for i in indices:
	ff = filenames[i - 1]	
	system("ulimit -t 500; ./" + master_cmd + " " + ff + "  > " + master_cmd + str(i) + ".res; mv master.out " + master_cmd + str(i) + ".out")
	system("echo " + str(grain_sizes[i]) + " `tail -1 " + master_cmd + str(i) + ".out | cut -f3-` >> master.dat")

# run pairs
system("rm -f pairs.dat")
for i in indices:
	ff = filenames[i - 1]
	system("ulimit -t 10000; ./pairs " + ff + " > pair" + str(i) + ".res; mv pairs.out pairs" + str(i) + ".out")
	system("echo " + str(grain_sizes[i]) + " `tail -1 pairs" + str(i) + ".out | cut -f3-` >> pairs.dat") 
