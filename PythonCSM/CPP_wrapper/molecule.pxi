cimport csmlib

cdef csmlib.python_molecule cppize_molecule(atoms, equivalence_classes):
    """ Convert the Python structures to a CPP-compatible molecule """
    cdef csmlib.python_molecule molecule
    cdef csmlib.python_atom bridge_atom

    for atom in atoms:
        bridge_atom.symbol = cs(atom.symbol)
        bridge_atom.adjacent = atom.adjacent
        bridge_atom.pos = atom.pos
        bridge_atom.mass = atom.mass
        molecule.atoms.push_back(bridge_atom)

    for ec in equivalence_classes:
        molecule.equivalenceClasses.push_back(list_to_vector_int(ec))

    return molecule

cdef pythonize_molecule(const csmlib.python_molecule &molecule):
    cdef int i;

    atoms = []
    for i in range(molecule.atoms.size()):
        atom = Atom(ps(molecule.atoms[i].symbol), vector_double_to_tuple(molecule.atoms[i].pos), useMass=False)
        atom.adjacent = vector_int_to_list(molecule.atoms[i].adjacent)
        atom._mass = molecule.atoms[i].mass

        atoms.append(atom)

    equivalenceClasses = []
    for i in range(molecule.equivalenceClasses.size()):
        oneClass = vector_int_to_list(molecule.equivalenceClasses[i])
        equivalenceClasses.append(oneClass)
    return (atoms, equivalenceClasses)
