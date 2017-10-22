#CSM Python API

##Introduction

A basic run of CSM would be as follows:
1. Read the molecule using the class `MoleculeReader`, which returns an instance of class 
`Molecule`
2. There are three calculation classes: 
`csm.calculations.Exact`, `csm.calculations.Trivial`, and `csm.calculations.Approx`.
Each has its own relevant arguments. creating an instance of the class automatically runs the
calculation in question, and the result can be accessed via the property `result`. 
In addition, certain statistics about the calculation will be available. Right now
the properties `dead_ends`, `perm_count`, and `num_branches` are available from Exact.
3. Additional information can be obtained from the CSMResult class that is returned with the property
`result`. For example, it is possible to call `compute_local_csm` and obtain the local CSM per atom
in the molecule
4. The result can be passed to a ResultWriter class, eg FileWriter

A simple example:

```
from csm.molecule.molecule import MoleculeReader
from csm.calculations import Exact, Operation
from csm.input_output.writers import FileWriter

mol=MoleculeReader('myfile.mol')
op=Operation("ch", sn_max=10)
calculation=Exact(op, mol)
FileWriter(calculation.result, "myfile.out", op_name, "mol", json_output=True)
```

This program will print a json dictionary of the results of an exact calculation on myfile.mol into the file myfile.out

##The Operation Class

this class is passed to all the calculation classes

###Methods

#### __init__(op, sn_max=8)

initalized with a case_insensitive op string, eg "C2", "c2", "s8", "S8", "CH", "ch", "CI", "cI", "Ci", "ci"

if sn_max is specified and the op is CH (or Ch, cH, ch), op order will be set to sn_max. 

When not specified, the default value is 8.

When not CH, it is ignored regardless.

##The MoleculeReader Class

A static class that creates instances of Molecule from files or strings.

###Methods:

####from_string:

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

####from_file:

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

####print_equivalence_class_summary

A method created for backwards compatibility. Initialization statistics regarding
equivalence classes that were previously printed during 
initialization, can be printed using this function.

It receives one argument, `display_chains`. When true, 
statistics on the molecule's chains will also be printed.

####to_json

returns a dictionary of the Molecule's field values, can be saved to a file and read from, 
using `from_json`

####from_json

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


##The ExactCalculation Class

todo

##The ApproxCalculation Class

todo

##The TrivialCalculation Class

todo

##The CSMResult Class

todo

##The ResultWriter Class

todo

##Additional Examples

For example, the `main()` of the main CSM commandline program is as follows:

```
    print("CSM version %s" % __version__)
    if not args:
        args = sys.argv[1:]
    csv_file = None
    #step one: parse args
    dictionary_args = get_split_arguments(args)
    try:
        #step two: read molecule from file
        mol=MoleculeReader.from_file(**dictionary_args)
        dictionary_args['molecule']=mol
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
            try:
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
            except TimeoutError:
                print("Timed out")
                return
        result=calc.result
        #step six: print the results
        if dictionary_args['calc_local']:
            result.compute_local_csm()
        r=FileWriter(result, **dictionary_args)
        r.write()
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