Guide to Permuters in CSM
=========================

Introduction
--------------

When calculating the exact CSM of a molecule, we go through each legal permutation of the molecules atoms,
and use this to calculate the direction of the reference plane, and the symmetry measure.

The molecule is divided into equivalence classes, atoms that are permitted to be permuted with each other.
Because it is not allowed to permute atoms between equivalence classes, each permutation within an 
equivalence class is independent of permutations within other equivalence classes.

This trait is taken advantage of in order to:
1. Recursively generate each possible legal permutation
2. Precalculate the values of the matrix A, vector B, and initial value of the CSM

The main permuter used is the *CythonPermuter*, which takes full advantage of #1, and is the main permuter 
used in all exact CSM calculations.

However, a secondary permuter, the *SinglePermPermuter*, is also used. The SinglePermPermuter does not 
generate any permutations. It receives a permutation, and uses the existing code supporting the CythonPermuter 
in order to precalculate the values of the matrix A, vector B, and initial value of the CSM, which are needed
further along in the claculation process.
The SinglePermPermuter is used in two cases:
1. The user has provided a permutation and wishes to measure the CSM of just that permutation
2. The approximate CSM algorithm, which generates estimated permutations and measures the CSM of each one.


The permuters use several component parts, and rather than beginning the explanation with the permuter, it
will be easier to understand if we begin with the components.


The Cache
-----------------

The CSM calculation makes repeated use of three vector operations on the atoms of the molecule:
the vector cross product, the vector inner product, and a sum of the outer product of (a,b) and (b,a)

In smaller molecules, it is efficient to store the results of these calculations rather than recalculating
them each time. For that purpose, we use the Cache class, which stores the results of each calculation for each
possible combination within each equivalence class.

A call to cache.inner_product(i,j) will retrieve from the Cache's stored inner_prpdouct values the result for
the atom indices (i,j)

In larger molecules, the Cache size quickly becomes so unwieldy as to be unusable. In addition, when calculating 
CSM for a single perm caching does not save any calculation time.
However, the code still makes calls to an object of type cache in order to get calculation results. 
To accomodate this, the FakeCache class was created.

The FakeCache is a class which inherits from Cache, but does not actually do any caching. For each call to 
FakeCache it calls the relevant calculation.

A call to fakecache.inner_product(i,j) will call the function inner_product on the vectors of atom i and atom j,
and return the result.


The PermCheckers
--------------------

Equivalence classes provide one constraint on legal permutations.

A second possible constraint is keepStructure. This allows only permutations which preserve the molecule's structure.

Example: 

     If atom A is bonded to atom B, and atom C is bonded to atom D.
     And (A,C) share an equivalence class and (B,D) share an equivalence class.

     keepStructure enforces the constraint the the permutation of A with C is only legal if B permutes with D,
     such that the bond between A and B is not broken.

     with the molecule: A-B=D-C

     C-D=B-A is a legal permutation
     C-B=D-A is not

within the CythonPermuter code that recursively builds permutations, each proposed added step of the permutation
is sent to a class PermChecker. 

There are two types of PermCheckers, both of which inherit from the "abstract" class PermChecker.

The *TruePermChecker* returns True without performing any check on the (to,from) atom indices sent to it.
The *StructurePermChecker* returns False if the proposed indices break the keepStructure constraint, and True otherwise.
It ascertains this by comparing between the neighbors of the (to) atom now and whether their permuted locations
share neighborship with the (from) atom. (if they have not yet been permuted, they are still potentially legal)

Which PermChecker is used is dependent on user-input arguments.


The CalcState
------------------

What is yielded by the permuter is not a permutation.
Rather, the permuter yields a CalcState, which contains a permutation, but also contains 
the values of the matrix A, the vector B, and the (not final) CSM.
These can then be used in the function calcRefPlane.
CalcState is a container class, which contains a convenience property, perm, for retrieving the permutation,
and a convenience function, copy, which copies the CalcState to another CalcState.


The PIPs (Perm in Progress)
---------------------

The PermutationInProgress is the base class PIP. 
It is used when we are only concerned with printing or counting permutations, 
and do not want to calculate anything else.

It contains the basic mechanism of building the perm in progress-- switch and unswitch.
Switch calles permchecker.is_legal and if so sets p[origin]=destination
For historical reasons, it also creates the reverse permutation and sets q[destination]=origin
Unswitch resets them back to -1.

In addition it has close_cycle and unclose_cycle, which are mainly used in its child class, PreCalcPIP. 
However, even in PermutationInProgress, close_cycle does copy the progress of the perm being built.

The main PIP, used for everything but counting/enumerating permutations, is PreCalcPIP.
It inherits unchanged switch and unswitch from its parent class.
However, it overwrites close_cycle and unclose_cycle:
on close cycle:
      PreCalcPIP calculates values for A, B, and CSM every time it closes a cycle, using the function partial_calculate.
      (each cycle's contributions to A,B,and CSM are independent from all other cycles)
	  In addition, it returns a copy of its saved calcstate before the calculation 
	  (which is saved by the parent function, and used in unclose_cycle)
on unclose cycle:
      The old_state is passed as an argument and saved as PreCalcPIP's calcState
	  
To faciliate the calculations, preCalcPIP also stores a small (size of the operation order) 
cache of the sines, cosines, and multipliers used in said calculation, 
which is created on initialization by the function _precalculate.
By default preCalcPIP uses a Cache, unless told otherwise (in which case it uses a FakeCache).


The Permuters
--------------------

The Permuter is the outwards facing class of the entire permutation generating enterprise.
It receives all the relevant arguments, and iteratively returns (via python's yield function)
for each permutation that has been calculated, a CalcState containing both the permutation and
that permutation's calculation values (A,B,CSM). 

#CythonPermuter

This has two outwards-facing functions.

The intialization function receives as arguments: 
     the molecule, the operation order, the operation type, 
     the argument keep_structure- if True, perm_checker= StructurePermChecker, otherwise perm_checker=TruePermChecker
	 and optional argument precalcute (default value true)- if false, will use PermInProgress instead of PreCalcPIP


The permute function receives no arguments.
It contains a pythonic yield loop, which calls on several inner functions.
_recursive_permute sends each group from the molecules equivalence groups to the function group_permuter
(and then recursively builds from the results of that group, and the remaining groups)
_group_permuter simply passes arguments to _group_recursive_permute (this does allow _recursive_permute to be more readable)

_group_recursive_permute is the main meat of the permutation algorithm. 
The algorithm is more easily explained within the code itself.


#SinglePermPermuter
The SinglePermPermuter contains within it a class called SinglePIP. 
SinglePIP inherits from PreCalcPIP, however, after initializing it calls close_cycle on the entire permutation.
This calculates A, B, and the CSM, and allows SinglePermPermuter to yield a single CalcState.
It is used whenever the A, B, and CSM of a signle permutation are desired- -in other words, --useperm and --approx
	 




