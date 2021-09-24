import numpy as np
from pymatgen.core import Molecule
from collections import Counter

def to_system_data(file_name) :
    mol = Molecule.from_file(file_name)
    elem_mol = list(str(site.species.elements[0]) for site in mol.sites)
    elem_counter = Counter(elem_mol)
    atom_names = elem_counter.keys()
    atom_numbs = elem_counter.values()
    atom_types = [list(atom_names).index(e) for e in elem_mol]
    
    system = {}
    system['atom_names'] = atom_names
    system['atom_numbs'] = atom_numbs
    system['atom_types'] = np.array(atom_types, dtype = int)
    center = mol.center_of_mass
    system['orig'] = np.array(center)

    mol.get_centered_molecule()
    system['coords'] = np.array([np.copy(mol.cart_coords)])

    return system
