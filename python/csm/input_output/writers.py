import json
from csm.input_output.formatters import format_CSM, non_negative_zero
import io
from openbabel import OBConversion
from csm.calculations.basic_calculations import check_perm_structure, check_perm_cycles
from csm.molecule.molecule import MoleculeReader


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

        # print dir

        f.write("\n DIRECTIONAL COSINES:\n\n")
        f.write("%lf %lf %lf\n" % (
            non_negative_zero(result.dir[0]), non_negative_zero(result.dir[1]), non_negative_zero(result.dir[2])))



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

        obmol = \
        MoleculeReader._obm_from_file(result.molecule._filename, result.molecule._format, result.molecule._babel_bond)[
            0]
        for to_remove in result.molecule._deleted_atom_indices:
            obmol.DeleteAtom(obmol.GetAtom(to_remove + 1))

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

        self.write_ob_molecule(obmol, format, f)

        # print resulting structure coordinates

        # update output coordinates to match symmetric structure
        for i in range(num_atoms):
            try:
                a = obmol.GetAtom(i + 1)
                a.SetVector(non_negative_zero(result.symmetric_structure[i][0]),
                            non_negative_zero(result.symmetric_structure[i][1]),
                            non_negative_zero(result.symmetric_structure[i][2]))
            except Exception as e:
                pass

        if format == 'pdb':
            f.write("\nMODEL 02")
        f.write("\nRESULTING STRUCTURE COORDINATES\n")
        self.write_ob_molecule(obmol, format, f)
        if format == 'pdb':
            f.write("END\n")

        # print dir

        f.write("\n DIRECTIONAL COSINES:\n\n")
        f.write("%lf %lf %lf\n" % (non_negative_zero(result.dir[0]),
                                   non_negative_zero(result.dir[1]),
                                   non_negative_zero(result.dir[2])))



    def write_ob_molecule(self, mol, format, f):
        """
        Write an Open Babel molecule to file
        :param mol: The molecule
        :param format: The output format
        :param f: The file to write output to
        :param f_name: The file's name (for extension-finding purpose)
        """
        conv = OBConversion()
        if not conv.SetOutFormat(format):
            raise ValueError("Error setting output format to " + format)

        # write to file

        try:
            s = conv.WriteString(mol)
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
    def __init__(self, result, op_name, format, print_local=False, json_output=False, *args, **kwargs):
        self.result = result
        self.op_name = op_name
        self.format = str.lower(format)
        self.print_local = print_local
        self.json_output = json_output
        self.result_string = self.get_result_string()

    def write(self):
        raise NotImplemented

    def get_result_string(self):
        result_io = io.StringIO()
        self._write_results(result_io)
        result_string = result_io.getvalue()
        result_io.close()
        return result_string

    def to_json(self):
        json_dict = {"Result":
            {
                "result_string": self.result_string,
                "molecule": self.result.molecule.to_json(),
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

        if self.print_local:
            self.write_local_csm(f)

        if self.op_name == "CHIRALITY":
            self.write_chirality(f)

        self.write_permutation(f)

    def write_header(self, f):
        f.write("%s: %.4lf\n" % (self.op_name, abs(self.result.csm)))
        f.write("SCALING FACTOR: %7lf\n" % non_negative_zero(self.result.d_min))

    def write_mol(self, f):
        if self.format == "csm":
            molwriter = CSMMolWriter()
        else:
            molwriter = OBMolWriter()

        molwriter.write(f, self.result, self.op_name, self.format)

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
            percent_structure = check_perm_structure(self.result.molecule, self.result.perm)
            print("The permutation found maintains",
                  str(round(percent_structure * 100, 2)) + "% of the original molecule's structure")

        except ValueError:
            print(
                "The input molecule does not have bond information and therefore conservation of structure cannot be measured")

        falsecount, num_invalid, cycle_counts = check_perm_cycles(self.result.perm, self.result.op_order, self.result.op_type)
        if True:  # falsecount > 0 or self.dictionary_args['calc_type'] == 'approx':
            print(
                "The permutation found contains %d invalid %s. %.2lf%% of the molecule's atoms are in legal cycles" % (
                    falsecount, "cycle" if falsecount == 1 else "cycles",
                    100 * (len(self.result.molecule) - num_invalid) / len(self.result.molecule)))
            for cycle_len in sorted(cycle_counts):
                valid = cycle_len == 1 or cycle_len == self.result.op_order or (
                    cycle_len == 2 and self.result.op_type == 'SN')
                count = cycle_counts[cycle_len]
                print("There %s %d %s %s of length %d" % (
                    "is" if count == 1 else "are", count, "invalid" if not valid else "",
                    "cycle" if count == 1 else "cycles",
                    cycle_len))

    def print_result(self):
        print("%s: %.6lf" % (self.op_name, abs(self.result.csm)))
        print("CSM by formula: %.6lf" % (self.result.formula_csm))

class FileWriter(ResultWriter):
    """
    A ResultWriter class that writes to a file 
    """
    def __init__(self, result, out_file_name, op_name, format, print_local=False, json_output=False, *args, **kwargs):
        self.out_file_name=out_file_name
        super().__init__(result, op_name, format(), print_local, json_output)
    def write(self):
        self.print_structure()
        self.print_result()
        if self.json_output:
            with open(self.out_file_name, 'w', encoding='utf-8') as f:
                json.dump(self.to_json(), f)
        else:
            with open(self.out_file_name, 'w', encoding='utf-8') as f:
                self._write_results(f)



class FolderWriter(ResultWriter):
    pass






def print_results(result, dictionary_args):
    """
    Prints the CSM calculation results
    :param result:
    :param dictionary_args:
    """
    filewriter = FileWriter(result, **dictionary_args)
    filewriter.write()
