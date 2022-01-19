README v0.2 / 28 March 2019

# CSM

## Introduction

The CSM (Continuous Symmetry Measure) program allows users to measure the continuous symmetry of molecules.
Molecules can be input in a variety of file formats (any supported by the OpenBabel program). The CSM can be measured
precisely for smaller molecules, or approximated for larger molecules where the runtime for calculating 
the exact measure would be unreasonable.

## Features

* The exact continuous symmetry measure calculator: The calculator goes through every single possible (valid) permutation
 of the molecule's atoms, and selects the permutation with the lowest continuous symmetry.
 
 * The approximate continuous symmetry measure calculator: For molecules large enough that going through all the permutations
 is not feasible, a good alternative is the approximate algorithm, which uses an iterative algorithm based on building
 permutations around a selection of possible symmetry axes to approximate the molecule's continous symmetry measure
 
 * The trivial continuous symmetry measure calculator: Particularly well suited to proteins with almost identical sub-polymers,
 the trivial continous symmetry measure returns the symmetry measure of the molecule without any displacement of molecules within
 each fragment.

## Installation

Installing csm cdurrently requires installing openbabel and openbabel python bindings first.
Contact the developers for more information.

## Usage

After installation, the program can be called from the command line. 

The `csm` program offers a choice of 5 commands: `read`, `exact`, `approx`,
`trivial`, `write`

the `read` and `write` commands enable CSM to be used as a pipe-able program.

`exact`, `approx`, and `trivial` allow use of the exact, approximate, and trivial CSM calculations.

A complete usage example would be as follows:

`csm exact c2 --input inputmolecule.mol --output outputfile.mol --keep-structure`

Help for any of the options can be accessed by entering `csm <OPTION_NAME> -h`, for example `csm exact -h`.

This will provide a full list of optional and required arguments, with explanations.

In addition to the possibility of using CSM from the command line, CSM can be accessed programmatically through its API, 
detailed in the file API.md (including examples)



## CITATIONS 

Please cite CSM using the following:

proteincsm:

> Tuvi-Arad I. and Alon G., "Improved Algorithms for Quantifying the Near Symmetry of Proteins: Complete Side Chains Analysis" (submitted). 
> Dryzun C., Zait A. and Avnir D., "Quantitative Symmetry and Chirality - A Fast Computational Algorithm for Large Structures: Proteins, Macromolecules, Nanotubes, and Unit Cells", J. Comput. Chem., 32, 2526-2538 (2011). 


Exact Algorithm for calculating the CSM of molecules:

> Alon G. and Tuvi-Arad I., "Improved algorithms for symmetry analysis: Structure preserving permutations", J. Math. Chem., 56(1), 193-212 (2018).

Original Code by Avnir and coworkers:

> Pinsky M., Dryzun C., Casanova D., Alemany P., Avnir D., "Analytical Methods for Calculating Continuous Symmetry Measures and the Chirality Measure", Journal of Computational Chemistry 29(16): 2712-2721 (2008).

> Zabrodsky H., Peleg S., Avnir D., "Continuous symmetry measures", Journal of the American Chemical Society 114(20): 7843-7851 (1992).



## CREDITS

**Science/Math:**

Dr. Inbal Tuvi-Arad, Department of Natural Sciences, The Open University of Israel

Dr. Gil Alon, Department of Mathematics and Computer Science, The Open University of Israel

Prof. David Avnir, Institute of Chemistry, The Hebrew University of Jerusalem

**Programming:**

The Research Software Company (researchsoftware.co.il)

**Testing, scripts and additional technical support:**

Sagiv Barhoom,The Open University of Israel

**Intensive testing:**

Yaffa Shalit, The Open University of Israel

The code for the hungarian algorithm is copyright (c) 2012, Jacob Frelinger


## CONTACT

For questions about the code, feature requests, and bug reports, feel free to use the csm google group: https://groups.google.com/g/csm-openu/


## License

This project is provided under the 3-clause BSD license:

Copyright (c) 2014, The Professor Avnir Group, The Hebrew University, Jerusalem
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY 
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.