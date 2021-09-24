import numpy as np
from pymatgen.core import Molecule
from collections import Counter

def to_system_data(file_name, protect_layer = 9) :
    mol = Molecule.from_file(file_name)
    elem_mol = list(str(site.species.elements[0]) for site in mol.sites)
    elem_counter = Counter(elem_mol)
    atom_names = list(elem_counter.keys())
    atom_numbs = list(elem_counter.values())
    atom_types = [list(atom_names).index(e) for e in elem_mol]
    natoms = np.sum(atom_numbs)
    print(atom_numbs)
    print(natoms)
    
    tmpcoord = np.copy(mol.cart_coords)
    cog = np.average(tmpcoord, axis = 0)
    dist = tmpcoord - np.tile(cog, [natoms, 1])
    max_dist = np.max(np.linalg.norm(dist, axis = 1))
    h_cell_size = max_dist + protect_layer
    cell_size = h_cell_size * 2
    shift = np.array([1,1,1]) * h_cell_size - cog
    tmpcoord = tmpcoord + np.tile(shift, [natoms, 1])
    tmpcell = cell_size * np.eye(3)

    system = {}
    system['atom_names'] = atom_names
    system['atom_numbs'] = atom_numbs
    system['atom_types'] = np.array(atom_types, dtype = int)
    center = mol.center_of_mass
    system['orig'] = np.array(center)

    mol.get_centered_molecule()
    system['coords'] = tmpcoord
    system['cells'] = tmpcell

    return system
