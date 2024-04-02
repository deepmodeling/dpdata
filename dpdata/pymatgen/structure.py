import numpy as np

try:
    from pymatgen.core import Structure
except ImportError:
    pass

def from_system_data(structure: Structure):
    symbols = [site.species_string for site in structure]
    atom_names = list(structure.symbol_set)
    atom_numbs = [symbols.count(symbol) for symbol in atom_names]
    atom_types = np.array([atom_names.index(symbol) for symbol in symbols]).astype(
            int
    )
    coords = structure.cart_coords
    cells = structure.lattice.matrix
    nopbc = not np.any(structure.pbc)

    info_dict = {
        "atom_names": atom_names,
        "atom_numbs": atom_numbs,
        "atom_types": atom_types,
        "coords": np.array([coords]),
        "cells": np.array([cells]),
        "nopbc": nopbc,
    }
    return info_dict