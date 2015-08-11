from argparse import ArgumentParser

from openbabel import OBAtomAtomIter, OBConversion, OBMol
from calculations.molecule import Atom, GetAtomicSymbol


def create_parser():
    parser = ArgumentParser()

    parser.add_argument('input', help='Input file')
    parser.add_argument('output', help='Output file')

    parser.add_argument('--babelbond', action='store_true', default=False, help='Let OpenBabel compute bonding')
    return parser


def read_file(parsed_args):
    conv = OBConversion()
    obmol = OBMol()
    ob_format = conv.FormatFromExt(parsed_args.input)
    if not ob_format:
        raise ValueError("Error discovering format from filename " + parsed_args.input)
    if not conv.SetInFormat(ob_format):
        raise ValueError("Error setting openbabel format")

    if not parsed_args.babelbond:
        conv.SetOptions("b", conv.INOPTIONS)
    if not conv.ReadFile(obmol, parsed_args.input):
        raise ValueError("Error reading file " + parsed_args.input + " using OpenBabel")
    return obmol


def write_output(obmol, output_file_name):
    """
    :param obmol: OBmol molecule
    :return: A list of Atoms
    """

    try:
        f = open(output_file_name, 'w', encoding='utf-8')
    except IOError:
        raise ValueError("Failed to open output file " + output_file_name + " for writing")

    size = obmol.NumAtoms()
    atoms = []

    for i in range(size):
        obatom = obmol.GetAtom(i + 1)
        # get symbol by atomic number
        symbol = GetAtomicSymbol(obatom.GetAtomicNum())
        position = (obatom.GetX(), obatom.GetY(), obatom.GetZ())

        atom = Atom(symbol, position, False)

        adjacent = []
        iter = OBAtomAtomIter(obatom)
        for neighbour_atom in iter:
            adjacent.append(neighbour_atom.GetIdx() - 1)
        atom.adjacent = adjacent

        atoms.append(atom)

    f.write("%i\n" % size)
    for i in range(size):
        f.write("%s%14lf %14lf %14lf\n" %
                (atoms[i].symbol,
                 atoms[i].pos[0],
                 atoms[i].pos[1],
                 atoms[i].pos[2]))

    for i in range(size):
        f.write("%d " % (i + 1))
        for j in atoms[i].adjacent:
            f.write("%d " % (j + 1))
        f.write("\n")


if __name__ == '__main__':
    # This code runs when you execute the script from PyCharm.
    parser = create_parser()
    parsed_args = parser.parse_args()  # Call the parse
    obmol = read_file(parsed_args)
    write_output(obmol, parsed_args.output)

