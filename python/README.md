The Python CSM Package
======================


Changes in version 0.12.1
------------------------
small changes to the arguments-- no limit on number for C and S anymore, fixed help slightly


Changes in version 0.12.0
------------------------
MAJOR CHANGE: moved from CythonPermuter (which calculated cycles by group) 
to ConstraintPermuter (which chooses the next atom to place based on constraints).

The new permuter is slower by around 2x than the old, because it is primarily written in not particularly efficient python.
However, because it traverses the permutation tree significantly more efficiently-- sometimes by several orders of magnitude-- it is 
capable of handling far more complex molecules and complex symmetries than the previous permuter (to say nothing of the original C++)

Minor changes/bugfixes:
1. Double bonds stored as only one bond (for calculating equivalence classes, etc)
2. Added the flag --no-constraint to use the original CythonPermuter instead of the new ConstraintPermuter
3. The flag --hungarian has been replaced with the flag --no-hungarian: in other words, hungarian is now default, and --no-hungarian indicates
that program should use the old, non-hungarian approximate algorithm
4. small fixes to printouts-- prints "CSM by formula" instead of "YAFFA CSM"
5. when adding bonds from pdb, _create_bondset was not called. This created an incorrect printout, but also meant bondset was inaccurate.




Changes in version 0.11.0
------------------------
MAJOR BUG FIX: correct symmetric structure is calculated at end of code (such that the program's CSM and the equation's CSM are the same)
minor aesthetic fixes: removed an unnecessary warning, added some scientific-notation printouts.


Changes in version 0.10.0
------------------------
1. New, faster, cpp-based Munkres
2. bugfix: printing separate models only for pdbs


Changes in version 0.9.3.1
------------------------
Added a printout of the result of yaffa_test

Changes in version 0.9.3
------------------------
features added:
1. added a function yaffa_test that calculates CSM from the final symmetric structure with somewhat different math from the CSM produced by the main algorithm
2. print normalization factor in normcsm
3. removed the flags --ignore-hy and --use-dir, at least for now
4. print the percentage of the original molecule's structure preserved in the permutation found
5. added an entire testing suite: loading from tests, running tests, etc
6. added an algorithm for choosing the next atom to use in building a cycle in the permutation
7. Add ENDMDL to models in pdb file output. (note: currently suspected to be buggy)
8. Added options to the approx calculation to reduce number of directions used: --use-best-dir, --no-orthogonal
9. Read connectivity from pdb files that have connectivity

bug fixes:
1. added chirality support to the approx calculation
2. redo create_Q on denormalize
3. bug in remove_hy where values of 0 were responding as if None
4. uncomment some printouts
5. various small cleanups





Changes in version 0.9.2
------------------------

MAJOR UPDATE

The input arguments (flags) have undergone a complete overhaul.
As a result, this version breaks backwards compatibility with other versions.
If you run this version with flags from the previous version it will not work in most cases (except where flags were incidentally left untouched, e.g. one-word flags)

The changes are as follows:
1. flags are now all in a single standardized format: entirely in lowercase, with words separated by dashes.
2. some defunct/non-supported flags have been removed. (no-limit, babel-test, time-only). they can be added back if we decide what we want them to do.
3. several changes have taken place regarding flags related to the approximate algorithm:
a. --findperm has been removed. --approx, which previously applied both --findperm and --findOutliers, is now equivalent to what --findperm was previously (in other words, --approx now calls the approximate algorithm, and does NOT also call findOutliers. use --find-outliers if you want to find outliers.)
b. --usedir is no longer a stand-alone argument: it is a modifier to --approx. if --usedir is called without approx, it will be ignored.
4. some calculation flags apply only to the exact calculation, others only to the approximate calculation. calling a non-applicable flag will print a warning to the user, but the program will continue with the specified type of calculation, ignoring the irrelevant flag.
Bugfixes:
fixed bug in --trivial that was causing error message and program exit
fixed bug with ignore-sym that caused it to not do anything
fixed bug in use-dir that caused it to not do anything
added in support for S1 symmetry (equivalent to CS)
fixed bug that caused sn-max argument to be ignored. when sn-max argument is not provided, default sn_max is 8.
fixed bug that caused error in --print-local



Changes in version 0.9.0
------------------------
Added the norm_csm commandline command
Added normalizations to norm_csm: 'standard', 'atom_number', 'fragment_mass_center', 'symmetric_fragment_mass_center'
All normalizations are currently untested, mathematically.

Changes in version 0.8.8
------------------------

modified printouts in print_approx: removed print of permutation, added comparison of old and new, tidied up
added a prinout describing the permutation tree (# of branches, # of dead ends)

Changes in version 0.8.7
------------------------
Numerical instability- order of operations in calculating matrices A and B led to differences in the sixteenth place after the decimal point.
These differences were compounded by a numerically unstable algorithm implemented in C++ for finding the roots of the sixth degree polynomial,
leading to differences in the 3rd or 4th place after the decimal.
The C++ implementation has been replaced with the slower, more stable numpy.roots, at the cost of around 15% program speed.

The symmetric structure is now calculated only at the very end of the approx algorithm (a minor improvement in code efficiency reflecting some code cleanup) 


Changes in version 0.8.6
------------------------
minor bug fixes to 0.8.5

Changes in version 0.8.5
------------------------

Features:
--print_approx flag has been added, to enable testing of Hungarian algorithm. when flag is used,
approx will print information about the direction, permutation, and csm of each iteration

if keepStructure or usechains are specified and no bonding information is provided, 
code will automatically use babelbond unless given flag --noBabel

c2 and s2 in approx now loop through chainperms 
(without --usechains, they, like all calculations, will simply calculate on a single chain)

Fixes:
polynomial coeffecients are now rounded to 13 places before being sent to calculation of root

int declarations were changed to long in order to fix compilation bug on openU servers

bug in adjustment of adjacency groups in remove_hy has been fixed

major bug in keep_structure's is_legal checker fixed.


Changes in version 0.8.4
------------------------
enabled --useperm with a permutation that swiches between atoms of 
different type or different connectivity. (prints message informing user)
fixed bugs in remove_hy/ignore_hy code
 

Changes in version 0.8.3
------------------------
Fix to bug in calculation of direction
-0 is printed as 0

Changes in version 0.8.2
------------------------
Added a slow Python implementation of the Hungarian Algorithm. Run with --hungarian

Changes in version 0.8.1
------------------------
Approx is working with properly with non-equivalent chains.
Approx is working again without chains.

Changes in version 0.8.0
------------------------
1. Approx is working


Changes in version 0.7.4
------------------------
1. Fixed the installation problems
2. Added the trivial calculation mode - for the trivial calculation.
   Run csm --trivial for trivial mode.

