import os
import re
from openbabel import OBConversion

from csm import __version__
from csm.calculations.basic_calculations import cart2sph
from csm.input_output.formatters import format_CSM, format_unknown_str, output_strings, non_negative_zero
from csm.molecule.molecule import MoleculeReader, MoleculeData, mol_string_from_obm


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

def get_line_header(index, result):
    index_str = "%02d" % (index + 1)  # start from 1 instead of 0
    return "L" + index_str + "_" + result.operation.op_code


class LegacyFormatWriter:
    def __init__(self, result, format):
        self.result = result
        self.format = format

    def write(self, f, write_local=False):
        f.write("%s: %s\n" % (self.result.operation.name, format_CSM(self.result.csm)))
        f.write("SCALING FACTOR: %7lf\n" % non_negative_zero(self.result.d_min))

        molecule_writer_class = OBMolWriter
        if self.format == "csm":
            molecule_writer_class = CSMFormatWriter
        molecule_writer = molecule_writer_class(self.result.molecule)
        if self.format == 'pdb':
            f.write("\nMODEL 01")
        f.write("\nINITIAL STRUCTURE COORDINATES\n")
        molecule_writer.write_original(f, self.result, self.format, legacy=True)

        if self.format == 'pdb':
            f.write("\nMODEL 02")
        f.write("\nRESULTING STRUCTURE COORDINATES\n")
        molecule_writer.write_symmetric(f, self.result, self.format, legacy=True)
        if self.format == 'pdb':
            f.write("END\n")

        f.write("\n DIRECTIONAL COSINES:\n\n")
        f.write("%lf %lf %lf\n" % (
            non_negative_zero(self.result.dir[0]), non_negative_zero(self.result.dir[1]),
            non_negative_zero(self.result.dir[2])))

        if write_local:
            self.write_local(f)

        if self.result.operation.name == "CHIRALITY":
            if self.result.op_type == 'CS':
                f.write("\n MINIMUM CHIRALITY WAS FOUND IN CS\n\n")
            else:
                f.write("\n MINIMUM CHIRALITY WAS FOUND IN S%d\n\n" % self.result.op_order)

        f.write("\n PERMUTATION:\n\n")
        for i in self.result.perm:
            f.write("%d " % (i + 1))
        f.write("\n")

    def write_local(self, f):
        sum = 0
        f.write("\nLocal CSM: \n")
        size = len(self.result.molecule.atoms)
        local_csm = self.result.local_csm
        for i in range(size):
            sum += local_csm[i]
            f.write("%s %7lf\n" % (self.result.molecule.atoms[i].symbol, non_negative_zero(local_csm[i])))
        f.write("\nsum: %7lf\n" % sum)

    def write_json(self, f):
        import io
        import json
        result_io = io.StringIO()
        self.write(result_io)
        result_string = result_io.getvalue()
        result_io.close()

        json_dict = {"Result":
            {
                "result_string": result_string,
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

        json.dump(json_dict, f)


class CSMFormatWriter:
    def __init__(self, molecule):
        self.molecule=molecule

    def write_symmetric(self, f, result, *args, **kwargs):
        self.write_molecule(f, result.symmetric_structure)

    def write_original(self, f, result, *args, **kwargs):
        self.write_molecule(f, result.molecule.atoms)

    def write_molecule(self, f, coordinates):
        size=len(coordinates)
        f.write("%i\n" % size)
        for i in range(size):
            f.write("%3s%10.5lf %10.5lf %10.5lf\n" %
                    (self.molecule.atoms[i].symbol,
                     non_negative_zero(coordinates[i][0]),
                     non_negative_zero(coordinates[i][1]),
                     non_negative_zero(coordinates[i][2])))

        for i in range(size):
            f.write("%d " % (i + 1))
            for j in self.molecule.atoms[i].adjacent:
                f.write("%d " % (j + 1))
            f.write("\n")

class OBMolWriter:
    def __init__(self, molecule):
        self.molecule=molecule
        self.metadata=molecule.metadata
        self.obmols=self.obm_from_molecule(molecule)
        self.obmol=self.obmols[0]
        if len(self.obmols)>1:
            raise ValueError("We are temporarily not supporting reading fragments from a file and writing results")
        self.data = MoleculeData(self.obmols[0])

    def obm_from_molecule(self, molecule):
        obmols = MoleculeReader._obm_from_strings(molecule.metadata.file_content,
                                                  molecule.metadata.format,
                                                  molecule.metadata.babel_bond)
        _atom_indices = []

        for mol_index, obmol in enumerate(obmols):
            num_atoms = obmol.NumAtoms()
            for atom_index in range(num_atoms):
                _atom_indices.append((mol_index, atom_index))

        for to_remove in molecule._deleted_atom_indices:
            mol_index, atom_index = _atom_indices[to_remove]
            obmol = obmols[mol_index]
            obmol.DeleteAtom(obmol.GetAtom(atom_index + 1))

        return obmols

    def set_obm_from_original(self, obmol, result):
        return self.set_coordinates(obmol, result.molecule.atoms)

    def set_obm_from_symmetric(self, obmol, result):
        return self.set_coordinates(obmol, result.symmetric_structure)

    def set_coordinates(self, obmol, coordinates):
        num_atoms = obmol.NumAtoms()
        for i in range(num_atoms):
            try:
                a = obmol.GetAtom(i + 1)
                a.SetVector(non_negative_zero(coordinates[i][0]),
                            non_negative_zero(coordinates[i][1]),
                            non_negative_zero(coordinates[i][2]))
            except Exception as e:
                print(e)
        return obmol

    def write_symmetric(self, f, result, format, legacy=False):
        obm=self.set_obm_from_symmetric(self.obmol, result)
        self.write_molecule(f, obm, format, legacy=legacy)

    def write_original(self, f, result, format, legacy=False):
        obm=self.set_obm_from_original(self.obmol, result)
        self.write_molecule(f, obm, format, legacy=legacy)

    def write_molecule(self, f, obmol, format, legacy=False):
        """
        Write an Open Babel molecule to file
        :param obmol: The molecule
        :param format: The output format
        :param f: The file to write output to
        :param legacy: replace end with endmdl
        """
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


    def replace_internal_index(self, string_to_modify):
        '''
        xyz files have an internal index unrelated to the ordinal index they have in file, this replaces them
        :param
        :return:
        '''
        start_index = string_to_modify.find("mol_index=")
        end_index = 0
        if start_index != -1:
            end_index = string_to_modify.find(";")
            start_index = start_index + 10
            inner_mol_index = string_to_modify[start_index:end_index]
            string_to_modify = string_to_modify.replace("mol_index=" + inner_mol_index, "mol_index=" + str(self.metadata.index+1))
        return string_to_modify

    def add_filename(self, string_to_modify):
        start_index = string_to_modify.find("mol_index=")
        end_index = 0
        if start_index != -1:
            end_index = string_to_modify.find(";")
        filename=self.metadata.filename
        if filename not in string_to_modify:
            title_start = string_to_modify[:end_index + 1]
            title_fin = string_to_modify[end_index + 1:]
            string_to_modify = title_start + "\tfilename=" + filename + title_fin
        return string_to_modify

    def modify_title(self, replace_internal_index=True, add_filename=True, add_is_symmetric=False):
        '''
        the title string, which can be used with obmol.SetTitle()
        :param replace_internal_index: replace the original indices in an xyz file with ordinal
        :param add_filename: add filename to header if not already present
        :param add_is_symmetric: add symmetric/original to title
        :return:
        '''
        title=self.metadata.initial_title
        if replace_internal_index:
            title = self.replace_internal_index(title)
        if add_filename:
            title = self.add_filename(title)
        if add_is_symmetric:
            pass
        return title


def mult_mol_writer(filename, format, result, index, symmetric=False, end=True):
        '''
        :param filename:
        :param obmols:
        :return:
        '''

        obmWriter=OBMolWriter(result.molecule)

        # get obmols
        obmols = obmWriter.obmols

        #requires result
        for i, obmol in enumerate(obmols):
            if symmetric:
                obmWriter.set_obm_from_symmetric(obmol, result)
            else:
                obmWriter.set_obm_from_original(obmol, result)

        # handle headers:
        metadata = result.molecule.metadata
        mol_name = metadata.appellation(no_file_format=True)
        mol_index = metadata.index + 1

        #requires result
        line_header = get_line_header(index, result)

        moleculewriters=[]
        for mol in obmols:
            title = obmWriter.modify_title()
            title += "\tSYM_TXT_CODE=" + line_header
            mol.SetTitle(title)

        # handle special cases
        if len(obmols) > 1 or not end:
            if format == "mol":
                string = ""
                for mol in obmols:
                    string += mol_string_from_obm(mol, format)
                    string += "\n$$$$\n"
                with open(filename, 'a') as file:
                    file.write(string)
                return

            elif format == "pdb":
                string = ""
                # It may be necessary to add "model        #" here (see: alternating)
                # but it's not clear, so I'm leaving it be until the topic comes up again
                for mol in obmols:
                    string += mol_string_from_obm(mol, format)
                string=string.replace("END", "ENDMDL")
                with open(filename, 'a') as file:
                    file.write(string)
                return

        # default, including for multiple obmols that aren't special case formats above
        with open(filename, 'a') as file:
            for mol in obmols:
                obmWriter.write_molecule(file, mol, format)

def alternating_mol_writer(results, filename, format):
        with open(filename, 'w') as file:
            file.write('')
        if format!="pdb":
            for mol_results in results:
                for index, result in enumerate(mol_results):
                    mult_mol_writer(filename, result, index, symmetric=False, end=False)
                    mult_mol_writer(filename, result, index, symmetric=True, end=False)
        else:
            i=0
            string=""
            for mol_results in results:
                for index, result in enumerate(mol_results):
                    obmWriter = OBMolWriter(result.molecule)
                    obmols = obmWriter.obmols

                    model_str = "MODEL        {}\n".format(i)
                    string+=model_str
                    for obmol in obmols:
                        obmWriter.set_obm_from_original(obmol, result)
                        for mol in obmols:
                            mol_string= mol_string_from_obm(mol, format)
                            modified = re.sub("MODEL        \d+", "", mol_string)
                            string+=modified
                    i += 1

                    model_str = "MODEL        {}\n".format(i)
                    string+=model_str
                    for obmol in obmols:
                        obmWriter.set_obm_from_symmetric(obmol, result)
                        for mol in obmols:
                            mol_string= mol_string_from_obm(mol, format)
                            modified=re.sub("MODEL        \d+", "", mol_string)
                            string+=modified
                    i += 1


            string = string.replace("END", "ENDMDL")
            string+="\nEND"
            with open(filename, 'w') as file:
                file.write(string)







class ScriptWriter:
    def __init__(self, results, format, out_file_name=None, polar=False, verbose=False, print_local=False, **kwargs):
        '''
        :param results: expects an array of arrays of CSMResults, with the internal arrays by command and the external
        by molecule. if you send a single CSM result or a single array of CSMResults, it will automatically wrap in arrays.
        a single array of results will be treated as multiple commands on one molecule
        :param format: molecule format to output to
        :param out_file_name: if none is provided, the current working directory/csm_results will be used
        '''
        try:
            if not isinstance(results[0], list):  # results is a single array
                results = [results]
        except TypeError:  # results isn't an array at all
            results = [[results]]

        self.verbose = verbose
        self.results = results
        self.format = format
        self.folder = out_file_name
        if not self.folder:
            self.folder = os.path.join(os.getcwd(), "csm_results")
        self.polar = polar
        self.print_local = print_local

    def format_CSM(self, result):
        if result.failed:
            return "n/a"
        return format_CSM(result.csm)

    def write(self):
        os.makedirs(self.folder, exist_ok=True)
        self.create_CSM_tsv(filename=os.path.join(self.folder, "csm.txt"))
        self.create_perm_tsv(filename=os.path.join(self.folder, "permutation.txt"))
        self.create_dir_tsv(filename=os.path.join(self.folder, "directional.txt"))
        # self.create_extra_tsv(filename = os.path.join(self.folder, "extra.tab"))
        self.create_extra_txt(filename = os.path.join(self.folder, "extra.txt"))
        if self.verbose:
            self.create_approx_statistics(out_folder = os.path.join(self.folder, 'approx'))
        self.create_legacy_files(out_folder = os.path.join(self.folder, 'old-csm-output'))
        self.create_symm_mols(filename = os.path.join(self.folder, "resulting_symmetric_coordinates." + self.format))
        self.create_initial_mols(filename = os.path.join(self.folder, "initial_normalized_coordinates." + self.format))
        self.write_version(filename = os.path.join(self.folder, "version.txt"))

    def write_version(self, filename):
        with open(filename, 'w') as file:
            file.write("CSM VERSION: " + str(__version__))

    def create_CSM_tsv(self, filename):
        # creates a tsv file with CSM per molecule
        with open(filename, 'w') as f:
            f.write("%-20s" % "#Molecule")
            for index, res in enumerate(self.results[0]):
                f.write("%-10s" % (get_line_header(index, res)))
            f.write("\n")
            for index, mol_results in enumerate(self.results):
                f.write("%-20s" % mol_results[0].molecule.metadata.appellation())
                for result in mol_results:
                    f.write("%-10s" % self.format_CSM(result))
                f.write("\n")

    def create_perm_tsv(self, filename):
        # creates a tsv for permutations (needs to handle extra long permutations somehow)
        with open(filename, 'w') as f:
            f.write("%-20s%-10s%-10s\n" % ("#Molecule", "#Command", "#Permutation"))
            for mol_index, mol_results in enumerate(self.results):
                for line_index, command_result in enumerate(mol_results):
                    f.write("%-20s" % (command_result.molecule.metadata.appellation()))
                    f.write("%-10s" % (get_line_header(line_index, command_result)))
                    write_array_to_file(f, command_result.perm, True)
                    f.write("\n")

    def create_dir_tsv(self, filename):
        with open(filename, 'w') as f:
            f.write("%-20s%-10s%10s%10s%10s\n" % ("#Molecule", "#Command", "X", "Y", "Z"))
            for mol_index, mol_results in enumerate(self.results):
                for line_index, command_result in enumerate(mol_results):
                    f.write("%-20s" % command_result.molecule.metadata.appellation())
                    f.write("%10s" % get_line_header(line_index, command_result))
                    write_array_to_file(f, command_result.dir, separator="%-10s")
                    f.write("\n")

    def create_extra_tsv(self, filename):
        # create headers
        # first row: cmd
        format_strings = []
        headers_arr_1 = []
        headers_arr_2 = []

        for line_index, command_result in enumerate(self.results[0]):
            format_string = ""
            for key in sorted(command_result.overall_statistics):
                format_string += "%-" + str(len(key) + 2) + "s"
                headers_arr_1.append(get_line_header(line_index, command_result))
                headers_arr_2.append(key)
            format_strings.append(format_string)
        format_strings[-1] += "\n"

        full_string = "".join(format_strings)

        with open(filename, 'w') as f:
            f.write("%-20s" % "#Molecule")
            f.write(full_string % tuple(headers_arr_1))
            f.write("%-10s" % " ")
            f.write(full_string % tuple(headers_arr_2))
            for mol_index, mol_results in enumerate(self.results):
                f.write("%-20s" % mol_results[0].molecule.metadata.appellation())
                for line_index, command_result in enumerate(mol_results):
                    f.write(
                        format_strings[line_index] % tuple([format_unknown_str(command_result.overall_statistics[key])
                                                            for key in sorted(command_result.overall_statistics)]))

    def create_approx_statistics(self, out_folder):
        # (in approx there's also a table with initial direction, initial CSM, final direction, final CSM, number iterations, run time, and stop reason for each direction in approx)
        # create headers
        # first row: cmd

        # first, check that at least one of the commands has running statistics:
        has_ongoing = False
        for command_result in self.results[0]:
            if command_result.ongoing_statistics:
                has_ongoing = True

        if not has_ongoing:
            return

        os.makedirs(out_folder, exist_ok=True)

        for mol_index, mol_results in enumerate(self.results):
            for line_index, command_result in enumerate(mol_results):
                name = command_result.molecule.metadata.appellation(no_file_format=True) + "_" + get_line_header(line_index,
                                                                                                                 command_result)
                filename = os.path.join(out_folder, name + ".tsv")
                if command_result.ongoing_statistics:
                    self._write_approx_statistics(filename, command_result.ongoing_statistics["approx"])

    def _write_approx_statistics(self, filename, stats):
        with open(filename, 'w') as f:
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

    def create_extra_txt(self, filename):
        # 2. the equivalence class and chain information (number and length)
        # 4. Cycle numbers and lengths
        with open(filename, 'w') as f:
            for item in output_strings:
                f.write(item)
                f.write("\n")

    def create_legacy_files(self, out_folder):
        os.makedirs(out_folder, exist_ok=True)
        for mol_index, mol_results in enumerate(self.results):
            for line_index, command_result in enumerate(mol_results):
                name = command_result.molecule.metadata.appellation(no_file_format=True) + "_" + get_line_header(line_index,
                                                                                                                 command_result)
                file_name = name + "." + command_result.molecule.metadata.format
                out_file_name = os.path.join(out_folder, file_name)
                try:
                    of = LegacyFormatWriter(command_result, self.format)
                    with open(out_file_name, 'w', encoding='utf-8') as f:
                        of.write(f)
                except Exception as e:
                    print("failed to write legacy file for" + file_name + ": " + str(e))

    def create_initial_mols(self, filename):
        # chained file of initial structures
        with open(filename, 'w') as file:
            file.write('')
        for mol_results in self.results:
            for index, result in enumerate(mol_results):
                if result.molecule.metadata.format == "csm":
                    writer = CSMFormatWriter(result.molecule)
                    with open(filename, 'a') as f:
                        writer.write_original(f, result)
                else:
                    mult_mol_writer(filename, self.format, result, index)

    def create_symm_mols(self, filename):
        # chained file of symmetric structures
        with open(filename, 'w') as file:
            file.write('')
        for mol_results in self.results:
            for index, result in enumerate(mol_results):
                if result.molecule.metadata.format == "csm":
                    writer = CSMFormatWriter(result.molecule)
                    with open(filename, 'a') as f:
                        writer.write_symmetric(f, result)
                else:
                    mult_mol_writer(filename, self.format, result, index, symmetric=True)

