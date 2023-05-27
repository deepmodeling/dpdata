import numpy as np

try:
    from pymatgen.core import Molecule
except ImportError:
    pass
from collections import Counter


def to_system_data(file_name, protect_layer=9):
    mol = Molecule.from_file(file_name)
    elem_mol = list(str(site.species.elements[0]) for site in mol.sites)
    elem_counter = Counter(elem_mol)
    atom_names = list(elem_counter.keys())
    atom_numbs = list(elem_counter.values())
    atom_types = [list(atom_names).index(e) for e in elem_mol]
    natoms = np.sum(atom_numbs)

    tmpcoord = np.copy(mol.cart_coords)

    system = {}
    system["atom_names"] = atom_names
    system["atom_numbs"] = atom_numbs
    system["atom_types"] = np.array(atom_types, dtype=int)
    # center = [c - h_cell_size for c in mol.center_of_mass]
    system["orig"] = np.array([0, 0, 0])

    system["coords"] = np.array([tmpcoord])
    system["cells"] = np.array([10.0 * np.eye(3)])
    return system
