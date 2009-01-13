#!/usr/bin/python

from os import system

st = 10 ** 0.5
grain_sizes = [10, 10 * st, 100, 100 * st, 1000, 1000 * st, 10000, 10000 * st, 100000]	
#grain_sizes = [100000 * st, 1000000, 1000000 * st, 10000000]
indices = [1,2,3,4,5]


filenames = []
f = open("sizes2" , "w")

i = 1
for gs in grain_sizes:
	f.write(str(i) + " " + str(gs) + "\n")
	system("./createNetworkForArticle.py " + str(gs) +"; mv params.net params" + str(i) + ".net" )
	filenames.append("params" + str(i) + ".net")
	i = i+1
f.close()

# run rate
for i in indices:
        ff = filenames[i - 1]
        system("ulimit -t 200; ./rate " + ff + " > rate" + str(i) + ".res; mv rate.out rate" + str(i) + ".out")

# run moment
for i in indices:
       ff = filenames[i - 1]
        system("ulimit -t 1000; ./moment " + ff + " > moment" + str(i) + ".res; mv moment.out moment" + str(i) + ".out")



# run pairs
for i in indices:
        ff = filenames[i - 1]
        system("ulimit -t 2000; ./pairs " + ff + " > pair" + str(i) + ".res; mv pairs.out pairs" + str(i) + ".out")

# run master
for i in indices:
        ff = filenames[i - 1]
        system("ulimit -t 1000000; ./master " + ff + " > master" + str(i) + ".res; mv master.out master" + str(i) + ".out")

