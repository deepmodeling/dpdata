import numpy as np
from pymatgen.core import Molecule
from collections import Counter
import dpdata

def to_system_data(file_name, protect_layer = 9) :
    mol = Molecule.from_file(file_name)
    elem_mol = list(str(site.species.elements[0]) for site in mol.sites)
    elem_counter = Counter(elem_mol)
    atom_names = list(elem_counter.keys())
    atom_numbs = list(elem_counter.values())
    atom_types = [list(atom_names).index(e) for e in elem_mol]
    natoms = np.sum(atom_numbs)
    
    tmpcoord = np.copy(mol.cart_coords)

    system = {}
    system['atom_names'] = atom_names
    system['atom_numbs'] = atom_numbs
    system['atom_types'] = np.array(atom_types, dtype = int)
    # center = [c - h_cell_size for c in mol.center_of_mass]
    system['orig'] = np.array([0, 0, 0])

    system['coords'] = [tmpcoord]
    system['cells'] = [10.0 * np.eye(3)]
    return system

def remove_pbc(system, protect_layer = 9):
    nframes = len(system["coords"])
    natoms = len(system['coords'][0])
    for ff in range(nframes):
        tmpcoord = system['coords'][ff]
        cog = np.average(tmpcoord, axis = 0)
        dist = tmpcoord - np.tile(cog, [natoms, 1])
        max_dist = np.max(np.linalg.norm(dist, axis = 1))
        h_cell_size = max_dist + protect_layer
        cell_size = h_cell_size * 2
        shift = np.array([1,1,1]) * h_cell_size - cog
        system['coords'][ff] = system['coords'][ff] + np.tile(shift, [natoms, 1])
        system['cells'][ff] = cell_size * np.eye(3)
    return system
