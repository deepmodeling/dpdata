from ase import Atom
from ase.db import connect
import numpy as np


def get_frames(fname, begin=0, step=1):
    asedb = connect(fname)
    num = asedb.count()
    at0 = asedb.get(1).toatoms()

    numbers = at0.numbers
    nat0 = at0[numbers.argsort()]
    chemical_symbols = nat0.get_chemical_symbols()
    unique_numbers = np.unique(numbers)
    unique_numbers.sort()
    atom_names =  [Atom(u).symbol for u in unique_numbers] 
    atom_numbs = [chemical_symbols.count(i) for i in atom_names]
    atom_types = np.array([atom_names.index(i) for i in chemical_symbols])

    all_coords = []
    all_cells = []
    all_energies = []
    all_forces = []
    all_virials = None

    for i in range(begin, num, step):
        ati = asedb.get(i+1)
        ats_ = ati.toatoms()
        numbers = ats_.numbers
        sorted_numbers = numbers.argsort()
        ats = ats_[sorted_numbers]
        data = ati.data
        energy = data['energy']
        cell = ats.get_cell().view()
        coord = ats.positions
        forces = data['forces'][sorted_numbers]
        all_coords.append(coord)
        all_cells.append(cell)
        all_forces.append(forces)
        all_energies.append(energy)

    return atom_names, atom_numbs, atom_types, np.array(all_cells), np.array(all_coords), np.array(all_energies), np.array(all_forces), all_virials
