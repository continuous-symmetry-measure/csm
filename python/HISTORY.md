The Python CSM Package
======================

Changes in version 1.0.0
------------------------
* parallel calculations with --parallel, with output during calculations
* check for cmd.txt before reading molecules
* copy cmd.txt to results folder
* --output-perms no longer receives a filename, it is a simple true/false flag that is saved to the results folder as perms.csv
* fixes to spacing in output
* cleaning up of chirality code for future easier debugging

Changes in version 0.22.4
------------------------
* bug fixes to writing csm files
* vast simplification of setup.py install process

Changes in version 0.22.3
------------------------
* primarily: adding a test suite
* clean up some of the arguments
* write initial molecule normalized
* bug fix: sn_max wasn't being used
* added filepath to metadata
* bug fixes to nromalization

Changes in version 0.22.2
------------------------
* fixed bug in pyx files
* initial_normalized_coordinates file are indeed normalized
* fixed bug in --print-denorm, and changed file output name
* fixed bug in --out-format
* slight tidying of arguments
* removed --print-branches
* replace permute-chains with dont-permute-chains
* added repr for CSMResult
* fixed bug in reading files with comfile from a folder
* fixed bug with --sn-max
* fixed outdated code in normalizations
* beautiful new test cases

Changes in version 0.22.1
------------------------
* bug fix select-atoms
* change default out folder
* change some error printouts
* remove old ScriptWriter, add WebWriter

Changes in version 0.22.0
------------------------
Major:
* printing results takes place as the calculation is done, rather than at the end
* added --parallel for parallezing across molecules, and changed existing parallel for approx to --parallel-dirs
Minor:
* renamed --legacy to --legacy-output, and added --legacy-files to toggle creation of legacy files in default output
* if more than 10 molecules are being calculated, result summaries are not printed to screen
* warning for missing connectivity is only printed for first molecule in a file containing many
* added help for comfile command
* increase molecule title space in output files
* change default result file when --output not specified, and change flag --not-unique to the clearer --overwrite
Bugs:
* additional bug fix in chirality

Changes in version 0.21.2
------------------------
* remove --print-local
* bug fix in chirality
* changes to molecule index printing
* add title with file information to molecule header

Changes in version 0.20.12
------------------------
Add option for alternating original/symm result file

Changes in version 0.20.9
------------------------

Bug fixes in output, approx, and some argument parsing

Changes in version 0.20.8
------------------------
* Compile-time check for Eigen version

Changes in version 0.20.7
-------------------------

* rename "command" to "comfile"
* expand parser help
* replace default pipe behavior with flag --pipe for exact, trivial, and approx
* move printing back to sys.stdout unless --pipe is specified
* save all printouts to a file
* added support for comments in command files, prefixed with #
* when a calculation fails on a specific molecule or command, continue, and report FailedResult for that line
* multiple tweaks to output, like changes in number of decimal places
* created class MoleculeMetaData for storing information about writing molecules unrelated to actual calculation.
includes adding support for molecule title, which can then be modified...
* handle select-mols with command. number according to original numbers, not new indices
* Many miscellaneous bug fixes, including symmetric structure in chirality not calculating properly
* add flag --not-unique for not creating new folders
* add flag --version for printing CSM version (usage: `csm --version`)
* add flag --verbose, without which approx spreadsheet is not created

Changes in version 0.20.5 (beta)
-------------------------

a variety of output changes are in the process of being added, but this version update 
is for the csm-website, and includes changes to the requirements

Changes in version 0.20.4 (beta)
-------------------------

* --connect flag added
* --no-overwrite flag will automatically generate a new name for a result file/folder, 
rather than overwriting existing
* --statistics flag removed
* default output and input directories added for --input, --input, --print-perms, command

I/O tweaks:
* enable support for molecule arguments in command lines
* output files renamed
* result tables in fixed-width
* save command.txt and version.txt 
* formatting for molecule name and command name changed
* resultwriter can accept single array or single result

Bug fixes from previous versions:

* csm format molecules now have their format saved internally
* csm format molecules now get written properly in results
* fixed bug with openbabel hanging forever on .txt files
* major bug causing all approx chirality calculations to be wrong fixed
* use-perm bug with operations fixed

Changes in version 0.20.1 (beta)
-------------------------
* --select-mols added
* --select-atoms added
* warn for invalid PDBs
* change to molecule json (and to result json, which contains a molecule json):
molecule no longer saves original file name, but rather saves original file content
* bug fixes


Changes in version 0.20.0 (beta)
-------------------------

**This update is partially backwards compatible**

* complete overhaul to output
* added reading multiple molecules from folder
* added command as an option (alongside exact, approx, trivial) which receives cmd.txt and 
runs and cms commands within, each command on its own line (eg "exact c2")
* --legacy preserves old output
* --simple prints only CSM
* for everything else, tables and files are created in a specified folder 
(if a file is provided, the folder will be extracted and used). inside specified folder,
legacy files and approx staitstics files are saved in their own folders
* bug fix to use-sequence
* suppress openbabel warnings


Changes in version 0.19.1
-------------------------
* added parallelization in approx calculation, called with --parallel and optional int argument for number of processes
* normcsm has received the same treatment as the main CSM program and had its arguments changed. 
* dircsm (test_directions) has been removed
* --normalize has been added to main CSM program, to allow calling normcsm from the main CSM program
* --statistics and --output-perms have default values for where the files are saved
* bug fix: reading from sys.stdin no longer causes program to hang
* various small bug fixes

Changes in version 0.19.0
-------------------------

**MAJOR UPDATE** This update is NOT backwards compatible.

The program has been split into parts, with commands `read`, `exact`, `approx`, `trivial`, and `write`.
This enables the program to be treated as separate, pipeable scripts
Arguments have been divided between their relevant commands,
although exact, approx, and trivial allow commands from read/write in conjunction with --input, --output

the flag --no-constraint has been removed


Changes in version 0.18.0
-------------------------
* switch to openbabel 2.4.1
* add flag --selective to approx
* add support for --statistics and --polar flag to the main csm program
* bug fix in normcsm

Changes in version 0.17.5
-------------------------
* fibonacci flag can receive  number of directions
* changes to printout from dircsm
* added saving of higher cycle maintenance to approx

Changes in version 0.17.4
-------------------------
Bug fixes

Changes in version 0.17.3
-------------------------
1. Removed random factor from fibonacci sphere
2. Added support for --keep-structure with --approx using the new permuter, DistanceConstraintPermuter. 
This Permuter sorts all the equivalence class placements by distance when rotated, and chooses constraints based on lowest distance
3. Added printouts explaining why we stopped iterating in --approx
4. Some speed improvements to permuter
5. internal code cleanup including bug fixes and fixes to outdated tests, upgrade to cython, etc
6. Can print statistics about each direction in dircsm using flag --statistics

Changes in version 0.17.2
-------------------------
* Added fibonacci sphere to dircsm
* removed argument --k from dircsm. instead there is --num-dirs, which can be used for either random-k or fibonacci-sphere
* added printout of chain permutation


Changes in version 0.17.1
-------------------------
Bug fixes

Changes in version 0.17.0
-------------------------
Creation of Code API:   
1. Read a molecule with class MoleculeReader.
Some information about the molecule's equivalence classes can be printed
via the function print_equivalence_class_summary()
2. Pass molecule and other arguments to 
calculations.Exact, calculations.Approx, or calculations.Trivial
3. Statistics from the calculation's run can be accessed, eg Exact.dead_ends
4. A CSMResult is returned. Further analysis like local_csm can be done on the results
5. The CSMResult and relevant arguments can be passed to a ResultWriter class, eg FileWriter

Added ability to save molecule to and load molecule from json, and json for results.

Removed the option to count perms

Added --remove-hy capabilities to --use-sequence

Stopping condition in approximate calculation has changed from 
proportional comparison of results to absolute comparison of difference

Changes in version 0.16.4
-------------------------
Added --timeout, receives number of seconds, with default of 5 minutes
If program runtime exceeds timeout, program exits

Bug fixes:
* fixed bug in dircsm that caused option 6 to crash
* fixed bug that caused Trivial calculation to crash
* fixed bug that caused Trivial to not handle ch correctly
* fixed formatting for prints in dircsm
* fixed bug that interpreted --k in dircsm as --keep-structure

Changes in version 0.16.3
-------------------------
* fixed bug when reading PDBs

Changes in version 0.16.2
-------------------------
* fix bug in use-sequence
* change direction choicesa from names to numbers: `0: user-input, 1: exact-structure, 2: greedy-first, 
 3: random-k, 4: cube-corners, 5: atom-vectors, 6: atom-vectors-orth`.
* added argument for dir_ouput file, which prints run output
* small changes to the printouts from test_dirs, ie print the best molecule isntead of the last 

Changes in version 0.16.1
-------------------------
* critical bug fix for --exact-structure and --greedy-first
* added flag --seed

Changes in version 0.16.0
-------------------------
 Added new script, `test_direction`. Can be run from command line with
 `run_direction <direction choice> symmetry in-file out-file additional-args`. 
 The available direction choices are `user-input, exact-structure, greedy-first, 
 random-k, cube-corners, atom-vectors, atom-vectors-orth`.
 Additional flags: `--dirs-file`,` --k`
 

Changes in version 0.15.5
-------------------------
* fixed bug that caused bondset not to be created/updated for pdbs
* fixed critical bug causing norm_csm to crash
* added priority for chain identification-- first alphabetic at 22, then numeric at 26, 
then fragment number if relevant
* fixed bug in remove-hy


Changes in version 0.15.4
-------------------------
* bug fixes to reading molecules


Changes in version 0.15.3
-------------------------
* added reading multiple molecules from a file ("fragments") as chains, with the argument --read-fragments
* added support for HETATM chain IDs when using --use-chains



Changes in version 0.15.2
-------------------------
* main change: normalization factors are now working as expected
* bug fixes like a bad print, program crashing on no chains.
* sightly improved printouts
* counting perms is done with constraints permuter unless --no-constraints is specified
* babel-bond no longer default. if there's no bonds with keep-structure, program exits. otherwise with no bonds program prints warning and continues


Changes in version 0.15.1.5
-------------------------
* reverse some changes to flags. --use-chains is back, --no-chains is gone. 
* got rid of the outdated flag --print-norm, now that we have an entire separate command for printing normalization information
* added the flag --output-norm FILE which prints various information from the normalization process to a file. may be removed in future versions.
* various changes to the normalization calculations, which are still buggy.

Changes in version 0.15.1
-------------------------
* changes to flags: --no-hungarian removed, classic algorithm is now accessed via --greedy flag, --new-chains is now --many-chains. --use-chains has been set to default, and replaced with the flag --no-chains to override the default instead. 
* babelbond is called unless --no-babel is specified, under all circumstances (previously only with --use-chains or --keep-structure)
* if normalization factors 1-4 are requested on a molecule without multiple fragments, a warning is printed to the user
* bug fix: normalization factor 5 has a square root taken from the factor

Changes in version 0.15.0
-------------------------
* Major change: Completely restructured internal structure of approx algorithm, to faciliate addition of new chains approx algorithm (and potentially other algorithms)
* added flag --new-chains. New chains refers to the new chains algorithm, which, given a direction, chooses the best chain permutation and then the best permutation using that chain permutation.
* fixed a bug in normcsm that was causing normalization option 3 to always fail
* some internal tidying up of the Molecule class, including adding/removing properties
* fixed a bug where chain_equivalences was always simply a list of all the chains
* added a series of tests to the testing suite that offer a more diverse set of options and molecule types than the large .xyz-based tests

Changes in version 0.14.1
------------------------
* added flag --json-output, for use in website, tests, etc. 
* if there was an error, it returns {"Error":"error message"} 
  otherwise it returns {"Result":{"csm":csm, "dir":dir, etc...}}

Changes in version 0.14.0
------------------------
* fixed several mathematical errors in normalization factors

Changes in version 0.13.9
------------------------
* Properly handle PDBs with TER lines when reading the PDB connectivity

Changes in version 0.13.8
------------------------
* When m_t_B_2 is (0,0,0), find the maximum lambda by looking at the lambdas, instead of finding the polynomial
 coefficients.

Changes in version 0.13.7
------------------------
* Display equivalence classes when using --use-sequence

Changes in version 0.13.6
------------------------
* Output csm molecule coordinates with 5 decimal places.
* Add the --output-branches flag
* Minor fixes to the --use-sequence flag

Changes in version 0.13.5
------------------------
* Compared floating point numbers a little less maticulously, fixing a bug on older compilers. 

Changes in version 0.13.4
------------------------
* Fix a bug that caused the exact calculation to skip some permutations.
* Improve the permuter speed, speeding up the exact calculation by a factor of 2.

Changes in version 0.13.2
------------------------
* Normalization factor flags have been changed to the numbers 0-6, in order to make typing them in less onerous. 
The help documentation under -h describes what normalization each number belongs to.

* Scientific notation printouts have been adjusted to 5 significant digits before the decimal point and 4 after.

* The --log flag has been removed.


Changes in version 0.13.1
------------------------

* Added flag --use-sequence for pdb molecules.
* Small bugfix: fixed bug causing --use-perm to crash when checking conservation of structure with no bond information


Changes in version 0.13.0
------------------------
* Normalization factors: 'standard', 'atom_number', 'fragment_center', 'symmetry_center', 'fragment_symm', 'fragment_perm', 'linear_csm'
have been added. (some of these existed in a preliminary form in the older version, however, there were many mathematical errors)

* The usage of the norm_csm program has changed.
instead of: norm_csm type input_molecule output_file normalization [additional arguments]
where normalization could only be a single type of normalization

  the new version is: norm_csm normalization type input_molecule output_file [additional arguments]
where as many normalizations as desired (of the available 7) can be specified.



Changes in version 0.12.1
------------------------
* Fixed bug in calls to perm from state that was causing 
singlepermpermuter and hence approx to crash

* Modified printout so molecule with no bonds was not reported as 
100% conserved

Changes in version 0.12.1
------------------------
* small changes to the arguments-- no limit on number for C and S anymore, fixed help slightly


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

