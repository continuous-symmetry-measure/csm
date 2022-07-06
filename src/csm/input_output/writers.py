import csv
import json
import os
import re
from pathlib import Path
import shutil

try:
    import openbabel.openbabel as ob
    from openbabel.openbabel import OBMol, OBConversion
except ImportError:
    import openbabel as ob
    from openbabel import OBMol, OBConversion



from csm import __version__
from csm.calculations.basic_calculations import cart2sph
from csm.input_output.formatters import format_CSM, format_unknown_str, output_strings, non_negative_zero
from csm.molecule.molecule import MoleculeReader




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
                f.write("%s" % ("%10.4lf" % item))


def get_line_header(index, operation):
    index_str = "%02d" % (index + 1)  # start from 1 instead of 0
    return "L" + index_str + "_" + operation.op_code


class LegacyFormatWriter:
    def __init__(self, result, format):
        self.result = result
        self.format = format

    # custom formatters so that even if changes are made elsewhere, legacy remains the same
    def format_CSM(self, csm):
        return "%.4lf" % (abs(csm))

    def non_negative_zero(self, number):
        import math
        if math.fabs(number) < 0.00001:
            return 0.0000
        else:
            return number

    def write(self, f, write_local=False):
        f.write("%s: %s\n" % (self.result.operation.name, self.format_CSM(self.result.csm)))
        f.write("SCALING FACTOR: %7lf\n" % self.non_negative_zero(self.result.d_min))

        molecule_writer = MoleculeWrapper(self.result,symmetric=False, normalized=True)
        if self.format == 'pdb':
            f.write("\nMODEL 01")
        f.write("\nINITIAL STRUCTURE COORDINATES\n")
        molecule_writer.clean_trait_titles()
        molecule_writer.write(f, self.result, consecutive=True)

        if self.format == 'pdb':
            f.write("\nMODEL 02")
        f.write("\nRESULTING STRUCTURE COORDINATES\n")
        molecule_writer.set_traits(symmetric=True, normalized=False)
        molecule_writer.clean_trait_titles()
        molecule_writer.write(f, self.result, consecutive=True)
        if self.format == 'pdb':
            f.write("END\n")

        f.write("\n DIRECTIONAL COSINES:\n\n")
        f.write("%lf %lf %lf\n" % (
            self.non_negative_zero(self.result.dir[0]), non_negative_zero(self.result.dir[1]),
            self.non_negative_zero(self.result.dir[2])))

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
            f.write("%s %7lf\n" % (self.result.molecule.atoms[i].symbol, self.non_negative_zero(local_csm[i])))
        f.write("\nsum: %7lf\n" % sum)



class MoleculeWriter:
    def __init__(self, molecule_wrapper):
        self.molecule_wrapper = molecule_wrapper
        self.format = molecule_wrapper.molecule.metadata.out_format
        self.write_molecule = self._write_obm_molecule
        if self.format == "csm":
            self.write_molecule = self._write_csm_molecule

    def write(self, f, coords, consecutive=False, model_number=None):
        self.write_molecule(f, coords, consecutive, model_number)

    def _write_csm_molecule(self, f, coordinates, *args, **kwargs):
        size = len(coordinates)
        f.write("%i\n" % size)
        for i in range(size):
            f.write("%3s%10.5lf %10.5lf %10.5lf\n" %
                    (self.molecule_wrapper.molecule.atoms[i].symbol,
                     non_negative_zero(coordinates[i][0]),
                     non_negative_zero(coordinates[i][1]),
                     non_negative_zero(coordinates[i][2])))

        for i in range(size):
            f.write("%d " % (i + 1))
            for j in self.molecule_wrapper.molecule.atoms[i].adjacent:
                f.write("%d " % (j + 1))
            f.write("\n")

    def _write_obm_molecule(self, f, coordinates, consecutive=False, model_number=None):
        """
        Write an Open Babel molecule to file
        :param obmol: The molecule
        :param format: The output format
        :param f: The file to write output to
        :param legacy: replace end with endmdl
        """
        if model_number is not None and self.format.lower() == "pdb":
            model_str = "MODEL     {}\n".format(model_number)
            f.write(model_str)

        obmols = self.molecule_wrapper.set_obm_coordinates(coordinates)

        if len(obmols) > 1:
            print("WARNING: result printing for fragments may have errors")

        for obmol in obmols:
            conv = OBConversion()
            if not conv.SetOutFormat(self.format):
                raise ValueError("Error setting output format to " + format)
            # write to file
            try:
                s = conv.WriteString(obmol)
            except (TypeError, ValueError, IOError):
                raise ValueError("Error writing data file using OpenBabel")

            if consecutive:
                if str.lower(self.format) == 'pdb':
                    s = re.sub(r"MODEL\s+\d+", "", s)
                    s = s.replace("END", "ENDMDL")
                if str.lower(self.format) in ['mol']:
                    s += "\n$$$$\n"
            f.write(s)


class MoleculeWrapper:
    class MoleculeData(object):
        """
        Taken from pybel https://github.com/openbabel/documentation/blob/master/pybel.py
        Store molecule data in a dictionary-type object

        Required parameters:
          `obmol` -- an Open Babel :obapi:`OBMol`
        Methods and accessor methods are like those of a dictionary except
        that the data is retrieved on-the-fly from the underlying :obapi:`OBMol`.
        Example:

        >>> mol = readfile("sdf", 'head.sdf').next() # Python 2
        >>> # mol = next(readfile("sdf", 'head.sdf')) # Python 3
        >>> data = mol.data
        >>> print data
        {'Comment': 'CORINA 2.61 0041  25.10.2001', 'NSC': '1'}
        >>> print len(data), data.keys(), data.has_key("NSC")
        2 ['Comment', 'NSC'] True
        >>> print data['Comment']
        CORINA 2.61 0041  25.10.2001
        >>> data['Comment'] = 'This is a new comment'
        >>> for k,v in data.items():
        ...    print k, "-->", v
        Comment --> This is a new comment
        NSC --> 1
        >>> del data['NSC']
        >>> print len(data), data.keys(), data.has_key("NSC")
        1 ['Comment'] False
        """

        def __init__(self, obmol):
            self._mol = obmol

        def _data(self):
            return [ob.toPairData(x) for x in self._mol.GetData() if
                    x.GetDataType() == ob.PairData or x.GetDataType() == ob.CommentData]

        def _testforkey(self, key):
            if not key in self:
                raise KeyError("'%s' is not a Key" % key)

        def keys(self):
            return [x.GetAttribute() for x in self._data()]

        def values(self):
            return [x.GetValue() for x in self._data()]

        def items(self):
            return zip(self.keys(), self.values())

        def __iter__(self):
            return iter(self.keys())

        def iteritems(self):
            return iter(self.items())

        def __len__(self):
            return len(self._data())

        def __contains__(self, key):
            return self._mol.HasData(key)

        def __delitem__(self, key):
            self._testforkey(key)
            self._mol.DeleteData(self._mol.GetData(key))

        def clear(self):
            for key in self:
                del self[key]

        def has_key(self, key):
            return key in self

        def update(self, dictionary):
            for k, v in dictionary.iteritems():
                self[k] = v

        def __getitem__(self, key):
            self._testforkey(key)
            return ob.toPairData(self._mol.GetData(key)).GetValue()

        def __setitem__(self, key, value):
            if key in self:
                pairdata = ob.toPairData(self._mol.GetData(key))
                pairdata.SetValue(str(value))
            else:
                pairdata = ob.OBPairData()
                pairdata.SetAttribute(key)
                pairdata.SetValue(str(value))
                self._mol.CloneData(pairdata)

        def __repr__(self):
            return dict(self.iteritems()).__repr__()

    def __init__(self, result, line_index=0, symmetric=False, normalized=False, original_title=""):
        self.result = result
        self.molecule = result.molecule
        self.line_index = line_index
        self.metadata = result.molecule.metadata
        self.format = self.metadata.format
        self._molecule_coords="uninitialized"
        if self.format != "csm":
            self.obmols = self.obms_from_molecule(self.molecule)
            self.obmol = self.obmols[0]
            if not symmetric:
                self.original_title_obmol = self.obmol.GetTitle()
            self.moleculedata = MoleculeWrapper.MoleculeData(self.obmol)
            self.set_initial_molecule_fields(original_title)
        else:
            self.original_title_obmol = None

        self.set_traits(symmetric, normalized)

    def write(self, file, model_number=0, consecutive=False):
        mw = MoleculeWriter(self)
        mw.write(file, self._molecule_coords, consecutive=consecutive, model_number=model_number)


    def set_traits(self, symmetric, normalized):
        self.symmetric=symmetric
        self.normalized=normalized
        self._molecule_coords=self.result.get_coords(symmetric, normalized)
        self._set_symmetric_title(symmetric)
        self._set_normalized_title(normalized)

    def obms_from_molecule(self, molecule):
        if self.format == "csm":
            return []
        if self.format == "pdb":
            return self.obms_from_pdb(molecule)
        if molecule.obmol:
            obmols = [molecule.obmol]
        else:
            obmols = MoleculeReader._obm_from_strings(molecule.metadata.file_content,
                                                  molecule.metadata.format,
                                                  molecule.metadata.babel_bond)
        
        self.build_obm_atom_indices(obmols)

        num_atoms_after_deleted = len(molecule.atoms)

        if len(self._obm_atom_indices) == num_atoms_after_deleted:  # the _deleted_atom_indices already deleted.
            return obmols
        for to_remove in reversed(molecule._deleted_atom_indices):
            mol_index, atom_index = self._obm_atom_indices[to_remove]
            obmol = obmols[mol_index]
            obmol.DeleteAtom(obmol.GetAtom(atom_index + 1))

        return obmols

    def build_obm_atom_indices(self, obmols):
        self._obm_atom_indices = []

        num_all_atoms_obmols = 0
        for mol_index, obmol in enumerate(obmols):
            num_atoms = obmol.NumAtoms()
            num_all_atoms_obmols += num_atoms
            for atom_index in range(num_atoms):
                self._obm_atom_indices.append((mol_index, atom_index))

    def obms_from_pdb(self, molecule):
        conv = OBConversion()
        conv.SetInFormat('pdb')
        conv.SetOptions("b", conv.INOPTIONS)
        obmols = []

        for content in molecule.metadata.file_content:
            lines = content.split('\n')
            fixed_string = ''

            for line in lines:
                tokens = line.split()
                if len(tokens) > 0 and tokens[0] in ['ATOM', 'HETATM']:
                    atom_index = int(tokens[1]) - 1 
                    if atom_index in molecule._deleted_atom_indices:
                        continue

                fixed_string += line + '\n'

            obmol = OBMol()
            conv.ReadString(obmol, fixed_string)
            obmols.append(obmol)

        self.build_obm_atom_indices(obmols)
        return obmols

    def set_obm_coordinates(self, coordinates):
        if self.format == "csm":
            return []
        num_atoms = len(self.molecule)
        for i in range(num_atoms):
            mol_index, atom_index = self._obm_atom_indices[i]
            obmol = self.obmols[mol_index]
            try:
                a = obmol.GetAtom(atom_index + 1)
                a.SetVector(non_negative_zero(coordinates[i][0]),
                            non_negative_zero(coordinates[i][1]),
                            non_negative_zero(coordinates[i][2]))
            except Exception as e:
                print(e)
        return self.obmols

    def insert_pdb_new_lines(self, string_to_modify):
        '''
        pdb lines cannot be longer than 78 characters, openbabel will automatically add new lines with correct header,
         if there's a newline
        :param string_to_modify:
        :return:
        '''
        if len(string_to_modify) < 77:
            return string_to_modify
        modified = ""
        curr_index = 0
        while len(string_to_modify[curr_index:]) > 77:
            string_to_modify = modified + string_to_modify[:77] + "\n" + string_to_modify[77:]
            curr_index = curr_index + 76
        return string_to_modify

    def append_description(self, description):
        if self.format == "csm":
            return
        if self.format in ["mol", "sdf"]:
            key = 'Comment'
            self.moleculedata[key] = description
        if self.format.lower() == "pdb":
            key = 'REMARK'
            description = "     " + description
            description = self.insert_pdb_new_lines(description)

            if key in self.moleculedata:
                description = self.moleculedata[key] + "\n" + description

            self.moleculedata[key] = description
        if self.format == "xyz":
            old_title = self.obmol.GetTitle()
            if description not in old_title:
                new_title = old_title + "  " + description
                self.obmol.SetTitle(new_title)

    def append_title(self, title):
        if self.format == "csm":
            return
        if self.format.lower() == "pdb":
            key = "TITLE"
            title = "     " + title
            title = self.insert_pdb_new_lines(title)

            if key in self.moleculedata:
                title = self.moleculedata[key] + "\n" + title

            self.moleculedata[key] = title
        else:
            old_title = self.obmol.GetTitle()
            new_title = old_title + "  " + title
            self.obmol.SetTitle(new_title)

    def set_title(self, title):
        if self.format == "csm":
            return
        if self.format.lower() == "pdb":
            key = "TITLE"
            title = "     " + title
            title = self.insert_pdb_new_lines(title)

            if key in self.moleculedata:
                title = self.moleculedata[key] + "\n" + title

            self.moleculedata[key] = title
        else:
            self.obmol.SetTitle(title)

    def clean_title(self, string_to_remove):
        if self.format == "csm":
            return
        original_title = self.obmol.GetTitle()
        if self.format.lower() == "pdb":
            original_title = self.moleculedata["TITLE"] if "TITLE" in self.moleculedata else ''

        if string_to_remove in original_title:
            new_title = original_title.replace(string_to_remove, "")
        else:
            return

        if self.format.lower() == "pdb":
            self.moleculedata["TITLE"] = new_title
        self.obmol.SetTitle(new_title)

    def set_initial_molecule_fields(self, original_title=""):
        if self.format == "csm":
            return
        original_title = original_title if original_title else self.obmol.GetTitle()
        new_title = original_title + "  "
        if self.metadata.appellation() not in original_title:
            new_title += self.metadata.appellation() + " "
        line_header = get_line_header(self.line_index, self.result.operation)
        if line_header not in original_title:
            new_title = new_title + line_header
        if new_title:
            self.set_title(new_title)
        description = "index=" + str(self.metadata.index + 1) + ";" + "filename: " + self.metadata.filename
        self.append_description(description)

    def _set_symmetric_title(self, symmetric=True):
        self.clean_title("(Original)")
        self.clean_title("(Symmetric)")
        if symmetric:
            self.append_title("(Symmetric)")
        else:
            self.append_title("(Original)")

    def _set_normalized_title(self, normalized=True):
        self.clean_title("Denormalized")
        self.clean_title("Normalized")
        if normalized:
            self.append_title("Normalized")
        else:
            pass #Denormalized is the default state and does not need to be particularly noted.
            #self.append_title("Denormalized")

    def clean_trait_titles(self):
        #used for legacy
        self.clean_title("(Original)")
        self.clean_title("(Symmetric)")
        self.clean_title("Denormalized")
        self.clean_title("Normalized")


def format_result_CSM(result):
    if result.failed:
        return "n/a"
    return format_CSM(result.csm)


class ContextWriter:
    def __init__(self, commands, out_format, out_file_name=None, *args, **kwargs):
        self.commands = commands
        self.out_format = out_format
        self.out_file_name = out_file_name

    def __enter__(self):
        return self

    def write(self, mol_result):
        raise NotImplementedError

    def __exit__(self, exc_type, exc_value, traceback):
        pass


class PipeContextWriter(ContextWriter):
    results_arr = []

    def write(self, mol_result):
        self.results_arr.append(mol_result)

    def __exit__(self, exc_type, exc_val, exc_tb):
        import sys
        sys.stdout.write(
            json.dumps([[result.to_dict() for result in mol_results_arr] for mol_results_arr in self.results_arr],
                       indent=4, sort_keys=True))


class LegacyContextWriter(ContextWriter):
    called = False

    def write(self, mol_result):
        if not self.out_file_name:
            raise ValueError("must provide ouput file for legacy writing")
        if len(mol_result) > 1 or self.called:
            raise ValueError("Legacy result writing only works for a single molecule and single command")
        writer = LegacyFormatWriter(mol_result[0], self.out_format)
        with open(self.out_file_name, 'w') as f:
            writer.write(f)
        self.called = True


class SimpleContextWriter(ContextWriter):
    def write(self, mol_result):
        for lin_index, line_result in enumerate(mol_result):
            print("mol", line_result.molecule.metadata.appellation(), "cmd", lin_index + 1, " CSM: ",
                  format_CSM(line_result.csm))


class ScriptContextWriter(ContextWriter):
    def __init__(self, commands, out_format,
                 out_file_name=None,
                 polar=False, verbose=False, print_local=False, argument_string="",
                 print_denorm=False,
                 legacy_files=False, json_output=False,
                 max_len_file_name=36, **kwargs):
        super().__init__(commands, out_format, out_file_name)
        self.molecule_index = 0  # used for writing pdb files
        self.result_index = 0  # used for writing pdb files
        self.verbose = verbose
        self.folder = out_file_name
        if not self.folder:
            self.folder = os.path.join(os.getcwd(), "csm_results")
        self.polar = polar
        self.print_local = print_local
        self.print_denorm = print_denorm
        self.argument_string = argument_string
        self.create_legacy_files = legacy_files
        self.create_json_file=json_output
        self.max_len_file_name = max_len_file_name + 4
        self.molecule_format = '%-' + str(self.max_len_file_name) + 's'
        self.original_title = None
        self.com_file=kwargs.get("command_file")
        self._kwargs=dict(kwargs)

        self.init_files()

    def get_line_header(self, index, operation):
        index_str = "%02d" % (index + 1)  # start from 1 instead of 0
        return "L" + index_str + "_" + operation.op_code

    def init_files(self):
        '''
        creates the initial files and writes their headers
        :return:
        '''
        os.makedirs(self.folder, exist_ok=True)

        self.csm_file = open(os.path.join(self.folder, "csm.txt"), 'w')
        self.csm_file.write(self.molecule_format % "#Molecule")
        for index, operation in enumerate(self.commands):
            self.csm_file.write("%10s" % (self.get_line_header(index, operation)))
        self.csm_file.write("\n")

        self.dir_file = open(os.path.join(self.folder, "directional.txt"), 'w')
        format_string = self.molecule_format + "%-10s%10s%10s%10s\n"
        self.dir_file.write(format_string % ("#Molecule", "#Command", "X", "Y", "Z"))

        self.perm_file = open(os.path.join(self.folder, "permutation.txt"), 'w')
        format_string = self.molecule_format + "%-10s%-10s\n"
        self.perm_file.write(format_string % ("#Molecule", "#Command", "#Permutation"))

        self.initial_molecules_file = open(
            os.path.join(self.folder, "initial_coordinates." + self.out_format),
            'w')
        self.symmetric_mols_file = open(os.path.join(self.folder, "resulting_symmetric_coordinates." + self.out_format),
                                        'w')

        # extra
        self.extra_file = open(os.path.join(self.folder, "extra.txt"), 'w')

        # approx
        if self.verbose:
            self.approx_folder = os.path.join(self.folder, 'approx')
            os.makedirs(self.approx_folder, exist_ok=True)
            self.trivial_folder = os.path.join(self.folder, 'trivial')
            os.makedirs(self.trivial_folder, exist_ok=True)

        # legacy
        if self.create_legacy_files:
            self.legacy_folder = os.path.join(self.folder, 'old-csm-output')
            os.makedirs(self.legacy_folder, exist_ok=True)

        #json
        if self.create_json_file:
            self.json_data=[]

        #version and commandline
        filename = os.path.join(self.folder, "version.txt")
        with open(filename, 'w') as file:
            if not self.argument_string.startswith('csm'):
                file.write('csm ' + self.argument_string)
            else:
                file.write(self.argument_string)
            file.write("CSM VERSION: " + str(__version__))

        #comfile
        if self.com_file:
            name=Path(self.com_file).name
            filename=os.path.join(self.folder, name)
            shutil.copy(self.com_file, filename)

    def create_perms_csv(self, args_dict, line_index):

        perms_csv_name = None
        output_perms = args_dict.get("output_perms", False)
        if output_perms:
            exact_folder = os.path.join(self.folder, 'exact')
            os.makedirs(exact_folder, exist_ok=True)
            filename = args_dict["molecule"].metadata.appellation(no_file_format=True) + "_" + get_line_header(
                        line_index, args_dict["operation"]) + ".csv"
            perms_csv_name = os.path.join(self.folder, 'exact', filename)
            csv_file = open(perms_csv_name, 'w')
            perm_writer = csv.writer(csv_file, lineterminator='\n')
            perm_writer.writerow(['op', 'Permutation', 'Direction', 'CSM'])
        return perms_csv_name

    def close_files(self):
        self.csm_file.close()
        self.dir_file.close()
        self.perm_file.close()
        self.initial_molecules_file.close()
        self.symmetric_mols_file.close()
        self.extra_file.close()

    def write_csm(self, mol_results):
        f = self.csm_file
        f.write(self.molecule_format % mol_results[0].molecule.metadata.appellation())
        for result in mol_results:
            f.write("%10s" % format_result_CSM(result))
        f.write("\n")

    def write_dir(self, mol_results):
        f = self.dir_file
        for line_index, command_result in enumerate(mol_results):
            f.write(self.molecule_format % command_result.molecule.metadata.appellation())
            f.write("%-10s" % get_line_header(line_index, command_result.operation))
            write_array_to_file(f, command_result.dir, separator="%-10s")
            f.write("\n")

    def write_perm(self, mol_results):
        f = self.perm_file
        for line_index, command_result in enumerate(mol_results):
            f.write(self.molecule_format % (command_result.molecule.metadata.appellation()))
            f.write("%-10s" % (get_line_header(line_index, command_result.operation)))
            write_array_to_file(f, command_result.perm, True)
            f.write("\n")

    def write_initial_mols(self, mol_results):
        molecule_wrapper = MoleculeWrapper(mol_results[0], symmetric=False, normalized=False)
        self.original_title = molecule_wrapper.original_title_obmol
        file = self.initial_molecules_file
        molecule_wrapper.write(file, consecutive=True, model_number=self.molecule_index + 1)


    def write_symmetric_mols(self, mol_results):
        file = self.symmetric_mols_file
        for command_index, command_result in enumerate(mol_results):
            if command_result.failed:
                continue
            molecule_wrapper = MoleculeWrapper(command_result, command_index, symmetric=True, normalized=False, original_title=self.original_title)
            molecule_wrapper.write(file, consecutive=True, model_number=self.result_index + 1)
            self.result_index += 1

    def write_legacy_files(self, mol_results):
        for line_index, command_result in enumerate(mol_results):
            name = command_result.molecule.metadata.appellation(no_file_format=True) + "_" + get_line_header(line_index,
                                                                                                             command_result.operation)
            file_name = name + "." + command_result.molecule.metadata.format
            out_file_name = os.path.join(self.legacy_folder, file_name)
            try:
                of = LegacyFormatWriter(command_result, self.out_format)
                with open(out_file_name, 'w', encoding='utf-8') as f:
                    of.write(f)
            except Exception as e:
                print("failed to write legacy file for" + file_name + ": " + str(e))

    def write_extra_txt(self, mol_results):
        # molecule.print_equivalence_class_summary(True)
        # result.print_summary(dictionary_args["legacy"])
        f = self.extra_file
        item = output_strings.fetch()
        while item is not None:
            f.write(item)
            f.write("\n")
            item = output_strings.fetch()


    def write_approx_file(self, mol_results):
        out_folder = self.approx_folder
        for line_index, command_result in enumerate(mol_results):
            name = command_result.molecule.metadata.appellation(no_file_format=True) + "_" + get_line_header(
                line_index,
                command_result.operation)
            filename = os.path.join(out_folder, name + ".tsv")
            if "approx" in command_result.ongoing_statistics:
                self._write_approx_statistics(filename, command_result.ongoing_statistics["approx"])

    def _write_approx_statistics(self, filename, stats_dict):
        with open(filename, 'w') as f:
            header_string="Op\tDir Index"
            if self.polar:
                header_string+="\tr_i\tth_i\tph_i"
            else:
                header_string+="\tx_i\ty_i\tz_i"
            header_string+="\tCSM_i"
            if self.polar:
                header_string+="\tr_f\tth_f\tph_f"
            else:
                header_string+="\tx_f\ty_f\tz_f"
            header_string+="\tCSM_f\tRuntime\t# Iter\tStop Reason"
            header_string+="\tChain Perm"
            header_string+="\tresult validity\tbest valid dir\tbest valid csm\tbest valid %"
            header_string+="\n"
            f.write(header_string)

            for op in stats_dict:
                stats=stats_dict[op]
                for direction_index, direction_dict in enumerate(stats):
                    dir = direction_dict['dir']
                    stat = direction_dict['stats']
                    start_str = op+"\t" + str(direction_index) + "\t"
                    try:
                        x, y, z = stat['start dir']
                        xf, yf, zf = stat['end dir']
                        if self.polar:
                            x, y, z = cart2sph(x, y, z)
                            xf, yf, zf = cart2sph(xf, yf, zf)
                        start_str = start_str + format_CSM(x) + "\t" + format_CSM(y) + "\t" + format_CSM(z) + "\t"

                        interim_str= start_str \
                                + format_CSM(stat['start csm']) + "\t" \
                                + format_CSM(xf) + "\t" + format_CSM(yf) + "\t" + format_CSM(zf) + "\t" \
                                + format_CSM(stat['end csm']) + "\t" \
                                + format_CSM(stat['run time']) + "\t" \
                                + str(stat['num iterations']) + "\t" \
                                + stat['stop reason'] + "\t"

                        interim_str+=stat["chain perm"]+"\t"
                        interim_str=interim_str+str(stat["validity"]["res valid"])+"\t"+str(stat["validity"]["dir"])+\
                                    "\t"+format_CSM(stat["validity"]["csm"])+"\t"+str(stat["validity"]["per"])
                        end_str=interim_str+"\n"
                        f.write(end_str)

                    except Exception as e:
                        try:
                            start_str = start_str + stat['stop reason'] + "\t"
                        finally:
                            f.write(start_str + "failed to read statistics:"+str(e)+"\n")

    def write_trivial_file(self, mol_results):
        out_folder = self.trivial_folder
        for line_index, command_result in enumerate(mol_results):
            name = command_result.molecule.metadata.appellation(no_file_format=True) + "_" + get_line_header(
                line_index,
                command_result.operation)
            filename = os.path.join(out_folder, name + ".tsv")
            if "trivial" in command_result.ongoing_statistics:
                stats=command_result.ongoing_statistics["trivial"]
                if "n/a" not in stats:
                    with open(filename, 'w') as f:
                        f.write("op\tchain perm\tcsm\tdir")
                        for op in stats:
                            for chain_perm in stats[op]:
                                chain_stats=stats[op][chain_perm]
                                f.write("\n"+op+"\t"+chain_perm+"\t"+format_CSM(chain_stats["csm"])+"\t"+str(chain_stats["dir"]))

    def write(self, molecule_results):
        # receives result array for single molecule, and appends to all the relevant files
        #print(self.folder)
        self.write_csm(molecule_results)
        self.write_dir(molecule_results)
        self.write_perm(molecule_results)
        self.write_initial_mols(molecule_results)
        self.write_symmetric_mols(molecule_results)
        self.write_extra_txt(molecule_results)
        if self.create_legacy_files:
            self.write_legacy_files(molecule_results)
        if self.verbose:
            self.write_approx_file(molecule_results)
            self.write_trivial_file(molecule_results)
        if self.create_json_file:
            self.json_data.append(molecule_results)
        self.molecule_index += 1

    def __exit__(self, exc_type, exc_value, traceback):
        '''
        Exit the runtime context related to this object.
        The parameters describe the exception that caused the context to be exited.
        If the context was exited without an exception, all three arguments will be None.
        If an exception is supplied, and the method wishes to suppress the exception (i.e., prevent it from being propagated), it should return a true value.
        Otherwise, the exception will be processed normally upon exit from this method.

        '''
        try:
            if self.out_format.lower() == "pdb":
                self.symmetric_mols_file.write("\nEND")
                self.initial_molecules_file.write("\nEND")
        except:  # no matter what, must close files
            pass

        self.close_files()

        try: #write json output
            if self.create_json_file:
                with open(os.path.join(self.folder,"json-results.json"), 'w') as f:
                    json.dump([[result.to_dict() for result in mol_results_arr] for mol_results_arr in self.json_data], f,
                               indent=4, sort_keys=True)
        except Exception as e: #not important enough to throw exception over
            print(e)


class WebWriter():
    '''
    receives results in an array of arrays, like the old ScriptWriter.
    retained for web-csm
    '''

    def __init__(self, results, format=None, out_folder=None, max_len_file_name=36, context_writer=ScriptContextWriter, **kwargs):
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
        self.results = results
        self.max_len_file_name = max_len_file_name + 4
        self.molecule_format = '%-' + str(self.max_len_file_name) + 's'


        self.commands = []
        for result in results[0]:
            self.commands.append(result.operation)

        if kwargs.get("out_format"):
            self.format = kwargs["out_format"]
        elif kwargs.get("in_format"):
            self.format = kwargs["in_format"]
        else:
            self.format = self.results[0][0].molecule.metadata.format
        self.folder = out_folder
        if not out_folder:
            if kwargs.get('out_file_name'):
                self.folder = kwargs['out_file_name']
            else:
                self.folder = os.path.join(os.getcwd(), "csm_results")
        self.kwargs = kwargs
        self.context_writer = context_writer

    def result_molecule_iterator(self):
        for mol_index, mol_results in enumerate(self.results):
            for command_index, command_results in enumerate(mol_results):
                yield MoleculeWrapper(command_results, mol_index, command_index)

    def write(self):
        self.writer = self.context_writer(self.commands, self.format, **self.kwargs)
        for mol_result in self.results:
            self.writer.write(mol_result)

    def create_extra_tsv(self, filename=None):
        if not filename:
            filename = os.path.join(self.folder, "extra.tab")
        # create headers
        # first row: cmd
        format_strings = []
        headers_arr_1 = []
        headers_arr_2 = []

        for line_index, command_result in enumerate(self.results[0]):
            format_string = ""
            for key in sorted(command_result.overall_statistics):
                format_string += "%-" + str(len(key) + 2) + "s"
                headers_arr_1.append(get_line_header(line_index, command_result.operation))
                headers_arr_2.append(key)
            format_strings.append(format_string)
        format_strings[-1] += "\n"

        full_string = "".join(format_strings)

        with open(filename, 'w') as f:
            f.write(self.molecule_format % "#Molecule")
            f.write(full_string % tuple(headers_arr_1))
            f.write("%-10s" % " ")
            f.write(full_string % tuple(headers_arr_2))
            for mol_index, mol_results in enumerate(self.results):
                f.write(self.molecule_format % mol_results[0].molecule.metadata.appellation())
                for line_index, command_result in enumerate(mol_results):
                    f.write(
                        format_strings[line_index] % tuple([format_unknown_str(command_result.overall_statistics[key])
                                                            for key in sorted(command_result.overall_statistics)]))

    def create_perm_tsv(self, filename=None):
        # preserved for csm-web
        if not filename:
            filename = os.path.join(self.folder, "permutation.txt")
        # creates a tsv for permutations (needs to handle extra long permutations somehow)
        with open(filename, 'w') as f:
            format_line = self.molecule_format + "%-10s%-10s\n"
            f.write(format_line % ("#Molecule", "#Command", "#Permutation"))
            for mol_index, mol_results in enumerate(self.results):
                for line_index, command_result in enumerate(mol_results):
                    f.write(self.molecule_format % (command_result.molecule.metadata.appellation()))
                    f.write("%-10s" % (get_line_header(line_index, command_result.operation)))
                    write_array_to_file(f, command_result.perm, True)
                    f.write("\n")

    def create_alternating_mols(self, filename=None):
        # preserved for csm web
        if not filename:
            filename = os.path.join(self.folder, "resulting_mols." + self.format)
        i = 1
        with open(filename, 'w') as file:
            for molecule_wrapper in self.result_molecule_iterator():
                molecule_wrapper.set_traits(symmetric=False, normalized=False)
                molecule_wrapper.write(file, consecutive=True, model_number=i)
                i += 1
                molecule_wrapper.set_traits(symmetric=True, normalized=False)
                molecule_wrapper.write(file, consecutive=True, model_number=i)
                i += 1
            if self.format.lower() == "pdb":
                file.write("\nEND")
