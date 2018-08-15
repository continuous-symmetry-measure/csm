import json

import os

from datetime import datetime

from csm import __version__
from csm.input_output.formatters import format_CSM, non_negative_zero, format_perm_count, format_unknown_str, \
    output_strings
import io
from openbabel import OBConversion
from csm.calculations.basic_calculations import check_perm_structure_preservation, check_perm_cycles, cart2sph
from csm.molecule.molecule import MoleculeReader, get_format, mol_string_from_obm
from csm.input_output.formatters import csm_log as print

from csv import DictWriter

def write_array_to_file(f, arr, add_one=False, separator=" "):
    '''
    :param f: filestream being written to
    :param arr: array being written
    :param add_one: some arrays should be written with 1-index instead of 0 index, if so add_one=True
    :param separator: default is " "
    '''
    for item in arr:


        if add_one:
            if item == "n/a":
                pass
            else:
                item = item + 1
            f.write(str(item) + separator)
        else:
            if item == "n/a":
                f.write("%10s" % item)
            else:
                f.write("%10s" % ("%.4lf" % item))

def print_structure(f, result):
            # print CSM, initial molecule, resulting structure and direction according to format specified
            try:
                percent_structure = check_perm_structure_preservation(result.molecule, result.perm)
                f.write("The permutation found maintains "+
                      str(round(percent_structure * 100, 2)) + "% of the original molecule's structure")

            except ValueError:
                f.write(
                    "The input molecule does not have bond information and therefore conservation of structure cannot be measured")

            if True:  # falsecount > 0 or self.dictionary_args['calc_type'] == 'approx':
                falsecount, num_invalid, cycle_counts, bad_indices = check_perm_cycles(result.perm, result.operation)
                f.write(
                    "The permutation found contains %d invalid %s. %.2lf%% of the molecule's atoms are in legal cycles" % (
                        falsecount, "cycle" if falsecount == 1 else "cycles",
                        100 * (len(result.molecule) - num_invalid) / len(result.molecule)))

                for cycle_len in sorted(cycle_counts):
                    valid = cycle_len == 1 or cycle_len == result.operation.order or (
                            cycle_len == 2 and result.operation.type == 'SN')
                    count = cycle_counts[cycle_len]
                    f.write("\nThere %s %d %s %s of length %d" % (
                        "is" if count == 1 else "are", count,
                        "invalid" if not valid else "",
                        "cycle" if count == 1 else "cycles",
                        cycle_len))

# molwriters
class CSMMolWriter:
    def write(self, f, result, op_name, format="csm"):
        self.print_output_csm(f, result, op_name)

    def print_output_csm(self, f, result, op_name):
        """
        Prints output in CSM format
        :param f: File to print to
        :param result: The result of the CSM calculation (a CSMState)
        :param calc_args: Calculation arguments to CSM
        """
        size = len(result.molecule.atoms)

        # print initial molecule

        f.write("\nINITIAL STRUCTURE COORDINATES\n%i\n\n" % size)
        for i in range(size):
            f.write("%3s%10.5lf %10.5lf %10.5lf\n" %
                    (result.molecule.atoms[i].symbol,
                     non_negative_zero(result.molecule.atoms[i].pos[0]),
                     non_negative_zero(result.molecule.atoms[i].pos[1]),
                     non_negative_zero(result.molecule.atoms[i].pos[2])))

        for i in range(size):
            f.write("%d " % (i + 1))
            for j in result.molecule.atoms[i].adjacent:
                f.write("%d " % (j + 1))
            f.write("\n")

        # print resulting structure coordinates

        f.write("\nMODEL 02 RESULTING STRUCTURE COORDINATES\n%i\n" % size)

        for i in range(size):
            f.write("%3s%10.5lf %10.5lf %10.5lf\n" %
                    (result.molecule.atoms[i].symbol,
                     non_negative_zero(result.symmetric_structure[i][0]),
                     non_negative_zero(result.symmetric_structure[i][1]),
                     non_negative_zero(result.symmetric_structure[i][2])))

        for i in range(size):
            f.write("%d " % (i + 1))
            for j in result.molecule.atoms[i].adjacent:
                f.write("%d " % (j + 1))
            f.write("\n")


def write_ob_molecule(obmol, format, f, legacy=False, header_append=""):
        """
        Write an Open Babel molecule to file
        :param obmol: The molecule
        :param format: The output format
        :param f: The file to write output to
        :param f_name: The file's name (for extension-finding purpose)
        """
        title=obmol.GetTitle()
        obmol.SetTitle(title+header_append)

        conv = OBConversion()
        if not conv.SetOutFormat(format):
            raise ValueError("Error setting output format to " + format)

        # write to file

        try:
            s = conv.WriteString(obmol)
        except (TypeError, ValueError, IOError):
            raise ValueError("Error writing data file using OpenBabel")

        if legacy:
            if str.lower(format) == 'pdb':
                s = s.replace("END", "ENDMDL")
        f.write(s)


class OBMolWriter:
    def write(self, f, result, op_name, format):
        self.legacy_print_output_ob(f, result, format, op_name)

    def legacy_print_output_ob(self, f, result, format, op_name):
        """
        Prints output using Open Babel
        :param f: File to write to
        :param result: The result of the CSM calculation (a CSMState)
        :param in_args: Input arguments to CSM
        :param calc_args: Calculation arguments to CSM
        :param out_args: Output arguments to CSM
        """
        # print initial molecule
        if format == 'pdb':
            f.write("\nMODEL 01")
        f.write("\nINITIAL STRUCTURE COORDINATES\n")

        obmols=self.obm_from_result(result)
        obmol=obmols[0]
        self.set_obm_from_original(obmol, result)
        write_ob_molecule(obmol, format, f, legacy=True)

        if format == 'pdb':
            f.write("\nMODEL 02")
        f.write("\nRESULTING STRUCTURE COORDINATES\n")

        self.set_obm_from_symmetric(obmol, result)
        write_ob_molecule(obmol, format, f, legacy=True)
        if format == 'pdb':
            f.write("END\n")

    def obm_from_result(self, result):
        obmols = MoleculeReader._obm_from_strings(result.molecule.metadata.file_content,
                                                  result.molecule.metadata.format,
                                                  result.molecule.metadata.babel_bond)
        self._atom_indices=[]

        for mol_index, obmol in enumerate(obmols):
            num_atoms=obmol.NumAtoms()
            for atom_index in range(num_atoms):
                self._atom_indices.append((mol_index, atom_index))

        for to_remove in result.molecule._deleted_atom_indices:
            mol_index, atom_index = self._atom_indices[to_remove]
            obmol=obmols[mol_index]
            obmol.DeleteAtom(obmol.GetAtom(atom_index + 1))

        return obmols

    def set_obm_from_original(self, obmol, result):
        num_atoms = obmol.NumAtoms()
        # update coordinates
        for i in range(num_atoms):
            try:
                atom = obmol.GetAtom(i + 1)
                atom.SetVector(non_negative_zero(result.molecule.atoms[i].pos[0]),
                               non_negative_zero(result.molecule.atoms[i].pos[1]),
                               non_negative_zero(result.molecule.atoms[i].pos[2]))
            except Exception as e:
                pass
        return obmol

    def set_obm_from_symmetric(self, obmol, result):
        num_atoms = obmol.NumAtoms()
        for i in range(num_atoms):
            try:
                a = obmol.GetAtom(i + 1)
                a.SetVector(non_negative_zero(result.symmetric_structure[i][0]),
                            non_negative_zero(result.symmetric_structure[i][1]),
                            non_negative_zero(result.symmetric_structure[i][2]))
            except Exception as e:
                pass
        return obmol




# resultwriters
class _ResultWriter:
    """
    A class for writing results. It can write molecules to various openbabel and CSM formats, 
    can write headers, local csm, permutation, etc. Most functions prefixed write_ expect a filestream. 
    The two print functions may be changed in the future but for now print to the screen
    Inheriting classes are recommended to inherit a write() function that can be called from main()
    """

    def __init__(self, result, format, print_local=False, *args, **kwargs):
        #make result an array of arrays if it not already so
        self.result = result
        self.op_name = result.operation.name
        self.format = str.lower(format)
        self.print_local = print_local
        self.result_string = self.get_result_string()

    def write(self):
        raise NotImplementedError

    def get_result_string(self):
        result_io = io.StringIO()
        self._write_results(result_io)
        result_string = result_io.getvalue()
        result_io.close()
        return result_string

    def to_dict(self):
        json_dict = {"Result":
            {
                "result_string": self.result_string,
                "molecule": self.result.molecule.to_dict(),
                "op_order": self.result.operation.order,
                "op_type": self.result.operation.type,
                "csm": self.result.csm,
                "perm": self.result.perm,
                "dir": list(self.result.dir),
                "d_min": self.result.d_min,
                "symmetric_structure": [list(i) for i in self.result.symmetric_structure],
                "local_csm": list(self.result.local_csm),
                "formula_csm": self.result.formula_csm,
                "normalized_molecule_coords": [list(i) for i in self.result.normalized_molecule_coords],
                "normalized_symmetric_structure": [list(i) for i in self.result.normalized_symmetric_structure],
            }
        }
        return json_dict

    def _write_results(self, f):
        self.write_header(f)
        self.write_mol(f)
        self.write_dir(f)

        if self.print_local:
            self.write_local_csm(f)

        if self.op_name == "CHIRALITY":
            self.write_chirality(f)

        self.write_permutation(f)

    def write_header(self, f):
        f.write("%s: %s\n" % (self.op_name, format_CSM(self.result.csm)))
        f.write("SCALING FACTOR: %7lf\n" % non_negative_zero(self.result.d_min))

    def write_mol(self, f):
        if self.format == "csm":
            molwriter = CSMMolWriter()
        else:
            molwriter = OBMolWriter()

        molwriter.write(f, self.result, self.op_name, self.format)

    def write_dir(self, f):
        f.write("\n DIRECTIONAL COSINES:\n\n")
        f.write("%lf %lf %lf\n" % (
            non_negative_zero(self.result.dir[0]), non_negative_zero(self.result.dir[1]),
            non_negative_zero(self.result.dir[2])))

    def write_local_csm(self, f):
        sum = 0
        f.write("\nLocal CSM: \n")
        size = len(self.result.molecule.atoms)
        for i in range(size):
            sum += self.result.local_csm[i]
            f.write("%s %7lf\n" % (self.result.molecule.atoms[i].symbol, non_negative_zero(self.result.local_csm[i])))
        f.write("\nsum: %7lf\n" % sum)

    def write_chirality(self, f):
        if self.result.op_type == 'CS':
            f.write("\n MINIMUM CHIRALITY WAS FOUND IN CS\n\n")
        else:
            f.write("\n MINIMUM CHIRALITY WAS FOUND IN S%d\n\n" % self.result.op_order)

    def write_permutation(self, f):
        f.write("\n PERMUTATION:\n\n")
        for i in self.result.perm:
            f.write("%d " % (i + 1))
        f.write("\n")

    def print_structure(self):
        try:
            percent_structure = check_perm_structure_preservation(self.result.molecule, self.result.perm)
            print("The permutation found maintains " +
                    str(round(percent_structure * 100, 2)) + "% of the original molecule's structure")

        except ValueError:
            print("The input molecule does not have bond information and therefore conservation of structure cannot be measured")

        falsecount, num_invalid, cycle_counts, bad_indices = check_perm_cycles(self.result.perm, self.result.operation)
        print(
            "The permutation found contains %d invalid %s. %.2lf%% of the molecule's atoms are in legal cycles" % (
                falsecount, "cycle" if falsecount == 1 else "cycles",
                    100 * (len(self.result.molecule) - num_invalid) / len(self.result.molecule)))

    def print_result(self):
        print("%s: %s" % (self.op_name, format_CSM(self.result.csm)))
        print("CSM by formula: %s" % (format_CSM(self.result.formula_csm)))

    def print_chain_perm(self):
        if len(self.result.molecule.chains)>1:
            print("Chain perm: ", self.result.chain_perm_string)


    def print_local(self, filename):
        result=self.result
        sum = 0
        with open(filename, 'w') as f:
            f.write("\nLocal CSM: \n")
            size = len(result.molecule.atoms)
            for i in range(size):
                sum += result.local_csm[i]
                f.write("%s %7lf\n" % (result.molecule.atoms[i].symbol, non_negative_zero(result.local_csm[i])))
            f.write("\nsum: %7lf\n" % sum)


class ApproxStatisticWriter:
    def __init__(self, statistics, stat_file_name, polar):
        self.statistics=statistics
        self.file_name = stat_file_name
        self.polar = polar


    def write(self):
        with open(self.file_name, 'w') as file:
            if self.polar:
                file.write("Index"
                        "\tr_i\tth_i\tph_i"
                       "\tCSM_i"
                       "\tr_f\tth_f\tph_f"
                       "\tCSM_f"
                        "\tRuntime"
                        "\t # Iter"
                        "\t Stop Reason"
                        "\n")
            else:
                file.write("Index\t"
                        "x_i\ty_i\tz_i"
                       "\tCSM_i"
                       "\tx_f\ty_f\tz_f"
                       "\tCSM_f"
                        "\tRuntime"
                        "\t # Iter"
                        "\t Stop Reason"
                        "\n")
            for index, key in enumerate(self.statistics):
                start_str=str(index) +"\t"
                try:
                    stat=self.statistics[key]
                    x,y,z=stat.start_dir
                    start_str =start_str + format_CSM(x)+ "\t"+format_CSM(y)+ "\t"+format_CSM(z)+ "\t"
                    xf,yf,zf=stat.end_dir
                    if self.polar:
                        x,y,z=cart2sph(x,y,z)
                        xf,yf,zf=cart2sph(xf,yf,zf)

                    file.write(start_str
                               + format_CSM(stat.start_csm)+"\t"
                               + format_CSM(xf) + "\t" + format_CSM(yf) + "\t" + format_CSM(zf) + "\t"
                               + format_CSM(stat.end_csm) + "\t"
                               + format_CSM(stat.run_time) + "\t"
                               + str(stat.num_iterations) + "\t"
                               + stat.stop_reason +
                               "\n")
                except:
                    file.write(start_str + "failed to read statistics\n")

class OldFormatFileWriter(_ResultWriter):
    """
    A _ResultWriter class that writes to a file
    """

    def __init__(self, result, out_file_name, print_local=False, json_output=False, out_format=None, *args, **kwargs):
        self.out_file_name = out_file_name
        self.json_output = json_output
        if not out_format:
            try:
                out_format = get_format(None, result.molecule.metadata.format)
            except ValueError:
                out_format = get_format(None, out_file_name)
        super().__init__(result, out_format, print_local)

    def write(self):
        self.print_structure()
        self.print_result()
        self.print_chain_perm()
        if self.json_output:
            with open(self.out_file_name, 'w', encoding='utf-8') as f:
                json.dump(self.to_dict(), f)
        else:
            with open(self.out_file_name, 'w', encoding='utf-8') as f:
                self._write_results(f)

def get_line_header(index, result):
    index_str = "%02d" % (index+1) #start from 1 instead of 0
    return "L"+index_str + "_" + result.operation.op_code

def get_mol_header(index, result):
    test= result.molecule.metadata.header()
    return test

class ScriptWriter:
    def __init__(self, results, format, out_file_name=None, polar=False, verbose=False, **kwargs):
        '''
        :param results: expects an array of arrays of CSMResults, with the internal arrays by command and the external
        by molecule. if you send a single CSM result or a single array of CSMResults, it will automatically wrap in arrays.
        a single array of results will be treated as multiple commands on one molecule
        :param format: molecule format to output to
        :param out_file_name: if none is provided, the current working directory/csm_results will be used
        '''
        try:
            if not isinstance(results[0], list): #results is a single array
                results = [results]
        except TypeError: #results isn't an array at all
            results=[[results]]

        self.verbose=verbose
        self.results=results
        self.format=format
        self.folder=out_file_name
        self.polar=polar

    def format_CSM(self, result):
        if result.failed:
            return "n/a"
        return format_CSM(result.csm)

    def write(self):
        os.makedirs(self.folder, exist_ok=True)
        self.create_CSM_tsv()
        self.create_perm_tsv()
        self.create_dir_tsv()
        #self.create_extra_tsv()
        self.create_extra_txt()
        if self.verbose:
            self.create_approx_statistics()
        self.create_legacy_files()
        self.create_symm_mols()
        self.create_initial_mols()
        self.write_version()

    def write_version(self):
        filename = os.path.join(self.folder, "version.txt")
        with open(filename, 'w') as file:
            file.write("CSM VERSION: "+str(__version__))

    def create_CSM_tsv(self):
    #creates a tsv file with CSM per molecule
        filename = os.path.join(self.folder, "csm.txt")
        with open(filename, 'w') as f:
            f.write("%-20s"%"#Molecule")
            for index, res in enumerate(self.results[0]):
                f.write("%-10s" % (get_line_header(index, res)))
            f.write("\n")
            for index, mol_results in enumerate(self.results):
                f.write("%-20s" % mol_results[0].molecule.metadata.header())
                for result in mol_results:
                    f.write("%-10s" % self.format_CSM(result))
                f.write("\n")

    def create_perm_tsv(self):
        #creates a tsv for permutations (needs to handle extra long permutations somehow)
        filename = os.path.join(self.folder, "permutation.txt")
        with open(filename, 'w') as f:
            f.write("%-20s%-10s%-10s\n" %("#Molecule", "#Command", "#Permutation"))
            for mol_index, mol_results in enumerate(self.results):
                for line_index, command_result in enumerate(mol_results):
                    f.write("%-20s" % (command_result.molecule.metadata.header()))
                    f.write("%-10s" % (get_line_header(line_index, command_result)))
                    write_array_to_file(f, command_result.perm, True)
                    f.write("\n")

    def create_dir_tsv(self):
        filename = os.path.join(self.folder, "directional.txt")
        with open(filename, 'w') as f:
            f.write("%-20s%-10s%10s%10s%10s\n" % ("#Molecule", "#Command", "X", "Y", "Z"))
            for mol_index, mol_results in enumerate(self.results):
                for line_index, command_result in enumerate(mol_results):
                    f.write("%-20s" % command_result.molecule.metadata.header())
                    f.write("%10s" % get_line_header(line_index, command_result))
                    write_array_to_file(f, command_result.dir, separator="%-10s")
                    f.write("\n")

    def create_extra_tsv(self):
        #create headers
        #first row: cmd
        format_strings=[]
        headers_arr_1=[]
        headers_arr_2=[]

        for line_index, command_result in enumerate(self.results[0]):
            format_string=""
            for key in sorted(command_result.overall_statistics):
                format_string+= "%-" + str(len(key)+2) + "s"
                headers_arr_1.append(get_line_header(line_index, command_result))
                headers_arr_2.append(key)
            format_strings.append(format_string)
        format_strings[-1]+="\n"

        full_string="".join(format_strings)

        filename = os.path.join(self.folder, "extra.tab")
        with open(filename, 'w') as f:
            f.write("%-20s" % "#Molecule")
            f.write(full_string % tuple(headers_arr_1))
            f.write("%-10s" % " ")
            f.write(full_string % tuple(headers_arr_2))
            for mol_index, mol_results in enumerate(self.results):
                f.write("%-20s" % mol_results[0].molecule.metadata.header())
                for line_index, command_result in enumerate(mol_results):
                    f.write(format_strings[line_index] % tuple([format_unknown_str(command_result.overall_statistics[key])
                                                                for key in sorted(command_result.overall_statistics)]))

    def create_approx_statistics(self):
        #(in approx there's also a table with initial direction, initial CSM, final direction, final CSM, number iterations, run time, and stop reason for each direction in approx)
        #create headers
        #first row: cmd

        #first, check that at least one of the commands has running statistics:
        has_ongoing=False
        for command_result in self.results[0]:
            if command_result.ongoing_statistics:
                has_ongoing=True

        if not has_ongoing:
            return

        out_folder=os.path.join(self.folder, 'approx')
        os.makedirs(out_folder, exist_ok=True)

        for mol_index, mol_results in enumerate(self.results):
            for line_index, command_result in enumerate(mol_results):
                name=command_result.molecule.metadata.header(no_file_format=True)+"_"+get_line_header(line_index, command_result)
                filename=os.path.join(out_folder, name+".tsv")
                if command_result.ongoing_statistics:
                    self._write_approx_statistics(filename, command_result.ongoing_statistics["approx"])

    def _write_approx_statistics(self, filename, stats):
        with open (filename, 'w') as f:
            if self.polar:
                f.write("Dir Index"
                        "\tr_i\tth_i\tph_i"
                        "\tCSM_i"
                        "\tr_f\tth_f\tph_f"
                        "\tCSM_f"
                        "\tRuntime"
                        "\t # Iter"
                        "\t Stop Reason"
                        "\n")
            else:
                f.write("Dir Index"
                        "\tx_i\ty_i\tz_i"
                        "\tCSM_i"
                        "\tx_f\ty_f\tz_f"
                        "\tCSM_f"
                        "\tRuntime"
                        "\t # Iter"
                        "\tStop Reason"
                        "\n")

            for direction_index, direction_dict in enumerate(stats):
                dir = direction_dict['dir']
                stat = direction_dict['stats']
                start_str = str(direction_index) + "\t"
                try:
                    x, y, z = stat['start dir']
                    start_str = start_str + format_CSM(x) + "\t" + format_CSM(y) + "\t" + format_CSM(z) + "\t"
                    xf, yf, zf = stat['end dir']
                    if self.polar:
                        x, y, z = cart2sph(x, y, z)
                        xf, yf, zf = cart2sph(xf, yf, zf)

                    f.write(start_str
                            + format_CSM(stat['start csm']) + "\t"
                            + format_CSM(xf) + "\t" + format_CSM(yf) + "\t" + format_CSM(zf) + "\t"
                            + format_CSM(stat['end csm']) + "\t"
                            + format_CSM(stat['run time']) + "\t"
                            + str(stat['num iterations']) + "\t"
                            + stat['stop reason'] + "\t"
                                                    "\n")
                except Exception as e:
                    try:
                        start_str = start_str + stat['stop reason'] + "\t"
                    finally:
                        f.write(start_str + "failed to read statistics\n")

    def create_extra_txt(self):
        #2. the equivalence class and chain information (number and length)
        #4. Cycle numbers and lengths
        filename = os.path.join(self.folder, "extra.txt")
        with open(filename, 'w') as f:
            for item in output_strings:
                f.write(item)
                f.write("\n")

    def create_legacy_files(self):
        out_folder=os.path.join(self.folder, 'old-csm-output')
        os.makedirs(out_folder, exist_ok=True)
        for mol_index, mol_results in enumerate(self.results):
            for line_index, command_result in enumerate(mol_results):
                name=command_result.molecule.metadata.header(no_file_format=True)+"_"+get_line_header(line_index, command_result)
                file_name=name+"."+command_result.molecule.metadata.format
                out_file_name= os.path.join(out_folder, file_name)
                try:
                    of=OldFormatFileWriter(command_result, out_file_name, out_format=command_result.molecule.metadata.format)
                    with open(out_file_name, 'w', encoding='utf-8') as f:
                        of._write_results(f)
                except Exception as e:
                    print("failed to write legacy file for"+file_name+": "+str(e))

    def ___create_legacy_file(self):
        #function for Devora to convert tests, not for lab
        out_file_name=os.path.join(self.folder, 'legacy.'+self.format)
        with open(out_file_name, 'a', encoding='utf-8') as f:
            for mol_index, mol_results in enumerate(self.results):
                for line_index, command_result in enumerate(mol_results):
                    of=OldFormatFileWriter(command_result, out_file_name, out_format=self.format)
                    of._write_results(f)

    def mult_mol_writer(self, filename, obmols, header=""):
        '''
        :param filename:
        :param obmols:
        :return:
        '''

        if len(obmols)>1:
            if self.format=="mol":
                string=""
                for mol in obmols:
                    string+=mol_string_from_obm(mol, self.format)
                    string+="\n$$$$\n"
                with open(filename, 'w') as file:
                    file.write(string)
                return

            elif self.format=="pdb":
                string=""
                for mol in obmols:
                    string+=mol_string_from_obm(mol, self.format)
                string.replace("END", "ENDMDL")
                string+="\nEND"
                with open(filename, 'w') as file:
                    file.write(string)
                return

        #default, including for multiple obmols that aren't special case formats above
        with open(filename, 'a') as file:
            for mol in obmols:
                write_ob_molecule(mol, self.format, file, header_append=header)


    def create_initial_mols(self):
        # chained file of initial structures
        filename = os.path.join(self.folder, "initial_normalized_coordinates." + self.format)
        with open(filename, 'w') as file:
            file.write('')
        for mol_results in self.results:
            for index, result in enumerate(mol_results):
                if result.molecule.metadata.format=="csm":
                    writer=CSMMolWriter()
                    with open(filename, 'a') as f:
                        writer.write(f, result, result.operation.name)
                    continue

                obmolwriter = OBMolWriter()
                obmols = obmolwriter.obm_from_result(result)
                for obmol in obmols:
                    obmolwriter.set_obm_from_original(obmol, result)
                header_append=" "+result.molecule.metadata.header()+"\tSYM_TXT_CODE="+get_line_header(index, result)
                self.mult_mol_writer(filename, obmols,header_append)

    def create_symm_mols(self):
        # chained file of symmetric structures
        filename = os.path.join(self.folder, "resulting_symmetric_coordinates." + self.format)
        with open(filename, 'w') as file:
            file.write('')
        for mol_results in self.results:
            for index, result in enumerate(mol_results):
                if result.molecule.metadata.format=="csm":
                    writer=CSMMolWriter()
                    with open(filename, 'a') as f:
                        writer.write(f, result, result.operation.name)
                    continue

                obmolwriter = OBMolWriter()
                obmols = obmolwriter.obm_from_result(result)

                for obmol in obmols:
                    obmolwriter.set_obm_from_symmetric(obmol, result)
                header_append=" "+result.molecule.metadata.header()+"\tSYM_TXT_CODE="+get_line_header(index, result)
                self.mult_mol_writer(filename, obmols,header_append)





