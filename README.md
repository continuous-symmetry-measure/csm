# CSM

# About

The CSM program calculates continuous symmetry and chirality measures of molecules with respect to a given point group G. Molecular coordinates can be provided in either xyz, sdf, mol or pdb format.
An online calculator is available at: https://csm.ouproj.org.il. 


## Features

* The codes cover the following point groups: Cs, Ci, Cn (n>=2), Sn (n= 2,4,6,…).
* Input structures can be in the form of a single molecule, concatenated file with multiple structures and a directory of molecules.
* When connectivity data is missing, OpenBabel is used to deduce it.

### Available commands
* comfile - Provide a command file for running calculations
* read - Read a molecule file into a json in CSM format
* write - Output the results of the calculation to a file. Must be used with piped input
* exact - Perform an exact CSM calculation for small-to-medium size molecules in xyz, mols, sdf and pdb file format. 
* approx - Approximate the CSM value. Relevant for protein homomers  in pdb file format. Partially tested for large molecules as well.
* trivial - Use the unit permutation to calculate the CSM for molecules and protein homomers.

## Citations

Please cite the CSM using the following:

### Exact algorithm:

> Alon G. and Tuvi-Arad I., "Improved algorithms for symmetry analysis: Structure preserving permutations", J. Math. Chem., 56(1), 193-212 (2018).

### Approx algorithm:

< Tuvi-Arad I. and Alon G., "Improved Algorithms for Quantifying the Near Symmetry of Proteins: Complete Side Chains Analysis", Journal of Cheminformatics, 11(1): 39 (2019).

< Dryzun C., Zait A. and Avnir D., "Quantitative Symmetry and Chirality - A Fast Computational Algorithm for Large Structures: Proteins, Macromolecules, Nanotubes, and Unit Cells", J. Comput. Chem., 32, 2526-2538 (2011).

### Original Code by Avnir and coworkers:

> Pinsky M., Dryzun C., Casanova D., Alemany P., Avnir D., "Analytical Methods for Calculating Continuous Symmetry Measures and the Chirality Measure", Journal of Computational Chemistry 29(16): 2712-2721 (2008).

<Zabrodsky, H.; Avnir, D. Continuous symmetry measures .4. Chirality. J. Am. Chem. Soc. 117: 462-473 (1995).

> Zabrodsky H., Peleg S., Avnir D., "Continuous symmetry measures", Journal of the American Chemical Society 114(20): 7843-7851 (1992).



## Usage

Input data requires a molecular geometry file and a choice of a point group
After installation, the program can be called from the command line. For example, to calculate the measure with respect to the C2 point group one should type:

'csm  exact c2 -- input input_mol.sdf --output output_dir --keep-structure'

For a list of all available options type `csm --help`

In addition to the possibility of using CSM from the command line, CSM can be accessed programmatically through its API. Details are in the file API.md (including examples).

##Installation

CSM can be installed on Windows and Linux machines.


### Build Instructions: Windows

Install [OpenBabel 3.1.1](https://github.com/openbabel/openbabel/releases/tag/openbabel-3-1-1)  
Test open babel with the command: `obabel -H` , if it doesn't work, try to restart your computer.  

From within the python folder, run:
`pip install -r requirements.txt`  

Run the build cython commands:  
`python\csm\CPP_wrapper> python .\setup.py build`  
`python\cython-munkres> .\rebuild.bat`  



### Build Instructions: Linux

???

Because installing openbabel correctly is a delicate and bug-prone process, an alternative method of installing CSM is available using Conda, and is described in the file conda_install_instructions.txt 

## Credits

**Science/Math:**

Prof. Inbal Tuvi-Arad, Department of Natural Sciences, The Open University of Israel

Dr. Gil Alon, Department of Mathematics and Computer Science, The Open University of Israel

Prof. David Avnir, Institute of Chemistry, The Hebrew University of Jerusalem

**Programming:**

The Research Software Company (researchsoftware.co.il)

**Testing, scripts and additional technical support:**

Sagiv Barhoom,The Open University of Israel

**Intensive testing:**

Yaffa Shalit, The Open University of Israel

The code for the hungarian algorithm is copyright (c) 2012, Jacob Frelinger


## Contact ##

For questions about the code, feature requests, and bug reports, feel free to use the CoSyM website users group at: https://groups.google.com/g/csm-openu. 

## License ##
This project is provided under the GPL 2 license. See `LICENSE.txt` for more information.
