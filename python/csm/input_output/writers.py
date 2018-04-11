import json

import os

from csm.input_output.formatters import format_CSM, non_negative_zero
import io
from openbabel import OBConversion
from csm.calculations.basic_calculations import check_perm_structure_preservation, check_perm_cycles, cart2sph
from csm.molecule.molecule import MoleculeReader, get_format
from csm.input_output.formatters import csm_log as print

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


class OBMolWriter:
    def write(self, f, result, op_name, format):
        self.print_output_ob(f, result, format, op_name)

    def print_output_ob(self, f, result, format, op_name):
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

        obmol=self.obm_from_result(result)
        obmol=self.set_obm_from_original(obmol, result)
        self.write_ob_molecule(obmol, format, f)

        if format == 'pdb':
            f.write("\nMODEL 02")
        f.write("\nRESULTING STRUCTURE COORDINATES\n")

        obmol=self.set_obm_from_symmetric(obmol, result)
        self.write_ob_molecule(obmol, format, f)
        if format == 'pdb':
            f.write("END\n")

    def obm_from_result(self, result):
        obmol = MoleculeReader._obm_from_file(result.molecule._filename,
                                              result.molecule._babel_bond,
                                              result.molecule._format)[0]
        for to_remove in result.molecule._deleted_atom_indices:
            obmol.DeleteAtom(obmol.GetAtom(to_remove + 1))
        return obmol

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



    def write_ob_molecule(self, obmol, format, f):
        """
        Write an Open Babel molecule to file
        :param obmol: The molecule
        :param format: The output format
        :param f: The file to write output to
        :param f_name: The file's name (for extension-finding purpose)
        """
        conv = OBConversion()
        if not conv.SetOutFormat(format):
            raise ValueError("Error setting output format to " + format)

        # write to file

        try:
            s = conv.WriteString(obmol)
        except (TypeError, ValueError, IOError):
            raise ValueError("Error writing data file using OpenBabel")
        if str.lower(format) == 'pdb':
            s = s.replace("END", "ENDMDL")
        f.write(s)


# resultwriters
class ResultWriter:
    """
    A class for writing results. It can write molecules to various openbabel and CSM formats, 
    can write headers, local csm, permutation, etc. Most functions prefixed write_ expect a filestream. 
    The two print functions may be changed in the future but for now print to the screen
    Inheriting classes are recommended to inherit a write() function that can be called from main()
    """

    def __init__(self, result, format, print_local=False, *args, **kwargs):
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
                "op_order": self.result.op_order,
                "op_type": self.result.op_type,
                "csm": self.result.csm,
                "perm": self.result.perm,
                "dir": list(self.result.dir),
                "d_min": self.result.d_min,
                "symmetric_structure": [list(i) for i in self.result.symmetric_structure],
                "local_csm": self.result.local_csm,
                "perm_count": self.result.perm_count,
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
        # print CSM, initial molecule, resulting structure and direction according to format specified
        try:
            percent_structure = check_perm_structure_preservation(self.result.molecule, self.result.perm)
            print("The permutation found maintains",
                  str(round(percent_structure * 100, 2)) + "% of the original molecule's structure")

        except ValueError:
            print(
                "The input molecule does not have bond information and therefore conservation of structure cannot be measured")

        if True:  # falsecount > 0 or self.dictionary_args['calc_type'] == 'approx':
            print(
                "The permutation found contains %d invalid %s. %.2lf%% of the molecule's atoms are in legal cycles" % (
                    self.result.falsecount, "cycle" if self.result.falsecount == 1 else "cycles",
                    100 * (len(self.result.molecule) - self.result.num_invalid) / len(self.result.molecule)))
            for cycle_len in sorted(self.result.cycle_counts):
                valid = cycle_len == 1 or cycle_len == self.result.op_order or (
                    cycle_len == 2 and self.result.op_type == 'SN')
                count = self.result.cycle_counts[cycle_len]
                print("There %s %d %s %s of length %d" % (
                    "is" if count == 1 else "are", count, "invalid" if not valid else "",
                    "cycle" if count == 1 else "cycles",
                    cycle_len))

    def print_result(self):
        print("%s: %s" % (self.op_name, format_CSM(self.result.csm)))
        print("CSM by formula: %s" % (format_CSM(self.result.formula_csm)))

    def print_chain_perm(self):
        if len(self.result.molecule.chains)>1:
            print("Chain perm: ", self.result.chain_perm_string())

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

class OldFormatFileWriter(ResultWriter):
    """
    A ResultWriter class that writes to a file 
    """

    def __init__(self, result, out_file_name, print_local=False, json_output=False, out_format=None, *args, **kwargs):
        self.out_file_name = out_file_name
        self.json_output = json_output
        if not out_format:
            try:
                out_format=get_format(None, out_file_name)
            except ValueError:
                out_format=get_format(None, result.molecule._filename)
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

class ScriptWriter:
    def __init__(self, results, format, commands, folder=os.getcwd()):
        self.results=results
        self.folder=folder
        self.commands=commands
        self.format=format

    def write(self):
        self.create_CSM_tsv()
        self.create_dir_tsv()
        self.create_initial_mols()
        self.create_output_txt()
        self.create_perm_tsv()
        self.create_result_out_folder()
        self.create_symm_mols()
        pass

    #logsymm1:
    def create_result_out_folder(self):
    # creates folder with each result.out for molecule
        pass

    def create_output_txt(self):
    #creates output.txt with all the outputs from screen
        pass

    def _write_initial_headers(self, f):
        f.write("\t")
        for index, command in enumerate(self.commands):
            f.write("L_"+str(index)+"\t")
            #args=command.split()
            #f.write("\t"+args[1])
        f.write("\n")

    def create_CSM_tsv(self):
    #creates a tsv file with CSM per molecule
        filename = os.path.join(self.folder, "csm.txt")
        with open(filename, 'w') as f:
            self._write_initial_headers(f)
            for mol_key, mol_results in self.results:
                f.write(mol_key)
                for result in mol_results:
                    f.write("\t"+format_CSM(result.csm))
                f.write("\n")

    #xyzsymm/pdbsymm:
    def create_dir_tsv(self):
        #creates a tsv for directions
        filename = os.path.join(self.folder, "dirs.txt")
        with open(filename, 'w') as f:
            self._write_initial_headers(f)
            for mol_key, mol_results in self.results:
                f.write(mol_key)
                for result in mol_results:
                    f.write("\t"+str(result.dir))
                f.write("\n")

    def create_perm_tsv(self):
        #creates a tsv for permutations (needs to handle extra long permutations somehow)
        filename = os.path.join(self.folder, "perms.txt")
        with open(filename, 'w') as f:
            self._write_initial_headers(f)
            for mol_key, mol_results in self.results:
                f.write(mol_key)
                for result in mol_results:
                    f.write("\t"+str(result.perm))
                f.write("\n")

    def create_initial_mols(self):
        # chained file of initial structures
        obmolwriter=OBMolWriter()
        filename=os.path.join(self.folder, "initial_normalized_coordinates."+self.format)
        with open(filename, 'w') as f:
            for mol_key, mol_results in self.results:
                for result in mol_results:
                    obmol=obmolwriter.obm_from_result(result)
                    obmol=obmolwriter.set_obm_from_original(obmol, result)
                    obmolwriter.write_ob_molecule(obmol, self.format, f)


    def create_symm_mols(self):
        # chained file of symmetric structures
        obmolwriter = OBMolWriter()
        filename = os.path.join(self.folder, "resulting_symmetric_coordinates." + self.format)
        with open(filename, 'w') as f:
            for mol_key, mol_results in self.results:
                for result in mol_results:
                    obmol = obmolwriter.obm_from_result(result)
                    obmol = obmolwriter.set_obm_from_symmetric(obmol, result)
                    obmolwriter.write_ob_molecule(obmol, self.format, f)




