README v0.0 / 15 March 2018

# CSM

## Introduction

The CSM (Continuous Symmetry Measure) program allows users to measure the continuous symmetry of molecules.
Molecules can be input in a variety of file formats (any supported by the OpenBabel program). The CSM can be measured
precisely for smaller molecules, or approximated for larger molecules where the runtime for calculating 
the exact measure would be unreasonable.

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
detailed in the file API.md

## Installation

CSM can be installed on Windows and Linux machines.

### Requirements
Before installing CSM, you must first install Openbabel (http://openbabel.org/wiki/Category:Installation), 
version 2.4.0 or later.

You must also `pip install both openbabel` (openbabel's python bindings) and `pip install numpy`

### Installation

CSM can be installed using:

`pip install csm --extra-index-url=https://repo.fury.io/theresearchsoftwarecompany`


## Credits

[to be added-- citation instructions]

Science/Math:

David Avnir (david.avnir@mail.huji.ac.il)

Inbal Tuvi-Arad (inbaltu@openu.ac.il)

Gil Alon (gil.alon@gmail.com)

Programming:

The Research Software Company (chelem.co.il)

Testing, scripts, additional technical support:

Sagiv

## Contact


### Contributing

You can contribute to the project by sending bug reports or feature requests to
Devora@chelem.co.il or via the contact form at chelem.co.il

### Help

Send requests for help to Devora@chelem.co.il or via the contact form at chelem.co.il

## License

This project is provided under the 3-clause BSD license.

The license should be in a separate file called LICENSE, so don't explain it in detail within your documentation. 

Also, don't forget to specify licenses of third-party libraries and programs you use.