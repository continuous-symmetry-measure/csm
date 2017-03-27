The files in this directory and sub-directories:

1. 1hw1.pdb and 1mwq.pdb are 2 dimers which the use-sequence option results a correct permutation and a low csm value.
However, running csm with the options --use-chains --trivial results a wrong permutation - the "real" trivial permutation (instead of "swiching the chains")- and a very high csm value.
The outputs can be found in the sub-dirs: use-sequence and use-chains-trivial

2. in the sub-dir length1 there are 2 dimers for which the csm value is low but there are circles of length 1. 
the run option was --use-sequence --approx.
we would like to understand what is wrong in the files
