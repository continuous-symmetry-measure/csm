# Changing cython-munkres to pypi munkres

## finding the args that use munkres
One of the pytest testing munkres: test_arguments.py\TestBasic\test_input_chain_perm

the args are:
approx c3 --input 1nc7.pdb  --use-sequence --use-chains --input-chain-perm CBDA2130.txt

## adding print to approx.pyx
this is where munkres is being used. line 245
print the matrix A before using munkres and then after using munkres, print the result.

## run csm with the args and save the output to a txt file
C:\Sources\CSM\csm\src\tests\argument_tests\files_for_tests [development] > ..\..\..\csm\main\csm_run.py approx c3 --input 1nc7.pdb  --use-sequence --use-chains --input-chain-perm CBDA2130.txt > C:\Sources\oldMunkres.txt

(I did this twice and compared the two output to make sure nothing is random or time affected.
Both outputs were identical.)

## Switch to PYPI munkres
install munkres from pypi: pip install munkres
edit approx.pyx to use the new munkres instead cython_munkres
rebuild
run csm again with the same args. This time, save the output to a different file:
C:\Sources\CSM\csm\src\tests\argument_tests\files_for_tests [development] > ..\..\..\csm\main\csm_run.py approx c3 --input 1nc7.pdb  --use-sequence --use-chains --input-chain-perm CBDA2130.txt > C:\Sources\newMunkres.txt

## Compare outputs
Compare oldMunkres and newMunkres using beyond compare.

## Result:
newMunkres is empty. But when output isn't being written to a file, there is some output. what is the reason to that?




