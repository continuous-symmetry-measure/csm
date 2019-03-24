#CSM Python API

##Introduction

A basic run of CSM would be as follows:
1. Read the molecule using the class `MoleculeReader`, which returns an instance of class 
`Molecule`
2. There are three calculation classes: 
`csm.calculations.Exact`, `csm.calculations.Trivial`, and `csm.calculations.Approx`.
Each has its own relevant arguments. All of them have the function calculate().
The result can be accessed via the property `result`, which contains an instance of the class `CSMResult`.
In addition, certain statistics about the calculation are available in the property `statistics`- the kind of
Statistics class returned varies by calculation.
3. Additional information can be obtained from the `CSMResult` class that is returned with the property
`result`. For example, it is possible to call `compute_local_csm` and obtain the local CSM per atom
in the molecule
4. The result can be passed to a `ResultWriter` class, eg `FileWriter`

A simple example:

```
from csm.molecule.molecule import MoleculeReader
from csm.calculations.data_classes import Operation
from csm.calculations import Exact
from csm.input_output.writers import FileWriter

mol=MoleculeReader('myfile.mol')
op=Operation("ch", sn_max=10)
calculation=Exact(op, mol)
calculation.calc()
FileWriter(calculation.result, "myfile.out", op, json_output=True)
```

This program will print a json dictionary of the results of an exact calculation on myfile.mol into the file myfile.out


##The MoleculeReader Class

A static class that creates instances of Molecule from files or strings.

###Methods:

####`from_string(...)`:

A static method that returns a Molecule instance read from a string

`from_string(string, format, initialize=True,
              use_chains=False, babel_bond=False,
             remove_hy=False, ignore_symm=False, use_mass=False):`

        :param string: the string to create the molecule from
        :param format: the format of the string (any BabelBond supported format, eg "mol", "xyz")
        :param initialize: boolean, default True, when True equivalence classes and chains are calculated for the molecule
        :param use_chains: boolean, default False, when True chains are read from the string
        :param babel_bond: boolean, default False, when True OpenBabel will attempt to guess connectivity information
        :param remove_hy: boolean, default False, when True hydrogen atoms will be removed
        :param ignore_symm: boolean, default False, when True all atom symbols will be "X" instead of elements
                            (this affects equivalence class calculation)
        :param use_mass: boolean, default False, when True atomic mass will be used, when False all atoms have a mass of 1
        :return: an instance of class Molecule

####`from_file(...)`:

A static method that returns a Molecule instance read from a file

`from_file(in_file_name, format=None, initialize=True,
                  use_chains=False, babel_bond=False,
                  remove_hy=False, ignore_symm=False, use_mass=False,
                  read_fragments=False, use_sequence=False,
                  keep_structure=False,
                  *args, **kwargs)`

```
        :param in_file_name: the name of the file to read the molecule from
        :param format: the format of the string (any BabelBond supported format, eg "mol", "xyz")
        :param initialize: boolean, default True, when True equivalence classes and chains are calculated for the molecule
        :param use_chains: boolean, default False, when True chains are read from the string
        :param babel_bond: boolean, default False, when True OpenBabel will attempt to guess connectivity information
        :param remove_hy: boolean, default False, when True hydrogen atoms will be removed
        :param ignore_symm: boolean, default False, when True all atom symbols will be "X" instead of elements
                            (this affects equivalence class calculation)
        :param use_mass: boolean, default False, when True atomic mass will be used, when False all atoms have a mass of 1
        :param read_fragments: boolean, default False, when True, multiple molecules in one file will be treated
                                as "chains" of one composite molecule
        :param use_sequence: boolean, default False, for pdb files only-- uses pdb sequence information to establish
                                equivalence (overwrites default equivalence)
        :param keep_structure: boolean, default False. If keep-structure is specified and molecule has no bond information,
                                the peogram will raise an error and exit
        :return: an instance of class Molecule
```

##The Molecule Class


Represents a molecule for CSM calculation.

###Methods:

####`print_equivalence_class_summary(display_chains=False)`

A method created for backwards compatibility. Initialization statistics regarding
equivalence classes that were previously printed during
initialization, can be printed using this function.

It receives one argument, `display_chains`. When true,
statistics on the molecule's chains will also be printed.

####`to_json()`

returns a dictionary of the Molecule's field values, can be saved to a file and read from,
using `from_json`

####`from_json(json)`

a static method, will return a Molecule with fields initialized
based on the values within a json dictionary, presumably created using
to_json.

###Properties:

####atoms
####norm_factor
####chains
####center_of_mass
####bondset
####equivalence_classes

##The Operation Class

This class is passed to all the calculation classes.

It is possible to replace this class with a NamedTuple. The only requirement for the object being passed is it
must have a field `type` and it must have a field `order`

###Methods

#### `__init__(op, sn_max=8)`

initalized with a case_insensitive op string, eg "C2", "c2", "s8", "S8", "CH", "ch", "CI", "cI", "Ci", "ci"

if sn_max is specified and the op is CH (or Ch, cH, ch), op order will be set to sn_max.

When not specified, the default value is 8.

When not CH, it is ignored regardless.

##The Calculation Class

ExactCalculation, ApproxCalculation, and TrivialCalculation all inherit from the base Calculation class.

All calculations expect to receive an operation and a molecule (and optionally additional arguments) when initializing.

All calculations have a function `calculate()` which performs the operation

All calculations have a property `result` which gives the CSMResult from the calculate() function
(is also generally returned directly by calls to that function )

ExactCalculation has a property `statistics` which returns an instance of `ExactStatistics`.

ApproxCalculation has a property `statistics` which returns an instance of `ApproxStatistics`.

TrivialCalculation has a filler property `statistics` which returns an empty dictionary.

##The ExactCalculation Class

###Initialization:

```
 def __init__(self, operation, molecule, keep_structure=False, perm=None, no_constraint=False, timeout=300, callback_func=None, *args, **kwargs):
        A class for running the exact CSM Algorithm
        :param operation: instance of Operation class or named tuple, with fields for name and order, that describes the symmetry
        :param molecule: instance of Molecule class on which the described symmetry calculation will be performed
        :param keep_structure: boolean, when True only permutations that maintain bond integrity will be measured 
        :param perm: a list of atom indiced describing one permutation of the molecule. Default None- When provided,
         only the provided permutation is measured
        :param no_constraint: boolean, default False, when False the constraints algorithm is used for the permuter, 
        when True the old permuter is used
        :param timeout: default 300, the number of seconds the function will run before timing out
        :param callback_func: default None, this function is called for every single permutation calculated with an argument of a single 
        CSMState, can be used for printing in-progress reports, outputting to an excel, etc.
```

###Additional properties:

####perm_count

The number of permutations yielded by the permuter

####dead_ends

The number of possible permutations that were cut off for breaking a constraint

####num_branches

Including dead-ends, the total number of decision branches in the permutation tree

##The ApproxCalculation Class
```
    def __init__(self, operation, molecule, approx_algorithm='hungarian', use_best_dir=False, get_orthogonal=True, detect_outliers=False, dirs=None, *args, **kwargs):
        """
        Initializes the ApproxCalculation
        :param operation: instance of Operation class or named tuple, with fields for name and order, that describes the symmetry
        :param molecule: instance of Molecule class on which the described symmetry calculation will be performed
        :param approx_algorithm: a string, can be 'hungarian', 'greedy', or 'many-chains'. Chooses the approximator alogrithm to use
        :param use_best_dir: use only the best direction 
        :param get_orthogonal: get orthogonal direction vectors from the main directions
        :param detect_outliers: detect outliers and use the imrpoved direction vectors
        :param dirs: a list of directions to use as initial dire (will override all other dir-choice options)
        """
```

###Additional methods

####The log() function

Throughout the approx calculation, calls are made to the ApproxCalculations log(*strings) function. 

The base class ApproxCalculation has a log function that does nothing, however, you can inherit from
ApproxCalculation and override the base behavior, for example printing the strings to the screen or
a file.


##The TrivialCalculation Class

```
    def __init__(self, operation, molecule, use_chains=True, *args, **kwargs):
        """
        :param operation: instance of Operation class or named tuple, with fields for name and order, that describes the symmetry.
        :param molecule: instance of Molecule class on which the described symmetry calculation will be performed.
        :param use_chains: default True. When True, all possible chain permutations with an identity perm on their components are measured.
                When false, only the pure identity perm is measured.
        """
```

##The CSMResult Class


This is a class returned from the program representing the results of the CSM calculation.



###properties

    """
    Properties:
    molecule: the original Molecule on which the calculation was called (denormalized)
    op_order: the order of the symmetry operation the calculation used
    op_type: the type of symmetry option (cs, cn, sn, ci, ch)
    csm: The final result CSM value
    perm: The final permutation that gave the result CSM value (a list of ints)
    dir: The final axis of symmetry
    perm_count: The number of permutations of the molecule
    local_csm: 
    d_min:
    symmetric_structure: the most symmetric structure when molecule permuted
    normalized_symmetric_structure: the normalized symmetric structure
    normalized_molecule_coords: the normalized molecule coordinates
    formula_csm: a recalculation of the symmetry value using a formula to measure symmetry. Should be close to identical
                    to algorithm result
    """

###methods

####create_symmetric_structure

automatically called in init, uses the calculated permutation and direction to create the nearest symmetric structure

####compute_local_csm

**Not** automatically called in init, when called, calculated the local csm (contribution of each atom to the overall CSM) 

##The ResultWriter Class

  A class for writing results. It can write molecules to various openbabel and CSM formats, 
    can write headers, local csm, permutation, etc. Most functions prefixed write_ expect a filestream. 
    The two print functions may be changed in the future but for now print to the screen
    
  The base class is ResultWriter, which generates the various things to be written to the provided stream.
  
  For now there is one inheriting class, FileWriter, which writes everything to a file (but eventually others will be added)
    
    
   The base args for ResultWriter are:
   
   `def __init__(self, result, op_name, format, print_local=False, *args, **kwargs):`
   
   For FileWriter:
   
   ` def __init__(self, result, out_file_name, op_name, format, print_local=False, json_output=False, *args, **kwargs):`
   
##Additional Examples

The `main()` of the main CSM commandline program is as follows:

```
        #step one: parse args
    dictionary_args = get_split_arguments(args)
    try:
        #step two: read molecule from file
        mol=MoleculeReader.from_file(**dictionary_args)
        dictionary_args['molecule']=mol
        dictionary_args['perm']=read_perm(**dictionary_args)
        #step three: print molecule printouts
        mol.print_equivalence_class_summary(dictionary_args['use_chains'])
        #step five: call the calculation
        if dictionary_args['calc_type'] == 'approx':
            if dictionary_args['print_approx']:
                class PrintApprox(Approx):
                    def log(self, *args, **kwargs):
                        print(*args)
                calc=PrintApprox(**dictionary_args)
            else:
                calc = Approx(**dictionary_args)
        elif dictionary_args['calc_type'] == 'trivial':
            calc = Trivial(**dictionary_args)
        else:
            # step four: create a callback function for the calculation
            # Outputing permutations
            csm_state_tracer_func= None
            if dictionary_args['perms_csv_name']:
                csv_file = open(dictionary_args['perms_csv_name'], 'w')
                perm_writer = csv.writer(csv_file, lineterminator='\n')
                perm_writer.writerow(['Permutation', 'Direction', 'CSM'])
                csm_state_tracer_func = lambda state: perm_writer.writerow(
                    [[p + 1 for p in state.perm],
                     state.dir,
                     state.csm, ])
            calc=Exact(**dictionary_args, callback_func=csm_state_tracer_func)
        try:
            calc.calc()
        except TimeoutError:
            print("Timed out")
            return
        result=calc.result
        #step six: print the results
        if dictionary_args['calc_local']:
            result.compute_local_csm()
        fw=FileWriter(result, **dictionary_args)
        fw.write()
        return result
    except Exception as e:
        if dictionary_args['json_output']:
            json_dict = {
                "Error": str(e)
            }
            with open(dictionary_args['out_file_name'], 'w', encoding='utf-8') as f:
                json.dump(json_dict, f)
        raise

    finally:
        if csv_file:
            csv_file.close()
```