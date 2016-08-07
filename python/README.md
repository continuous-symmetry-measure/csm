The Python CSM Package
======================
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

