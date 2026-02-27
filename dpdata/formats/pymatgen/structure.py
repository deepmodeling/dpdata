from __future__ import annotations

import numpy as np


def from_system_data(structure) -> dict:
    """Convert one pymatgen structure to dpdata's datadict."""
    symbols = [ii.specie.symbol for ii in structure]
    atom_names = list(structure.symbol_set)
    atom_numbs = [symbols.count(symbol) for symbol in atom_names]
    atom_types = np.array([atom_names.index(symbol) for symbol in symbols]).astype(int)
    coords = structure.cart_coords
    cells = structure.lattice.matrix
    if all(structure.pbc):
        pbc = True
    elif not any(structure.pbc):
        pbc = False
    else:
        raise ValueError(f"Partial pbc condition {structure.pbc} is not supported")

    info_dict = {
        "atom_names": atom_names,
        "atom_numbs": atom_numbs,
        "atom_types": atom_types,
        "coords": np.array([coords]),
        "cells": np.array([cells]),
        "orig": np.zeros(3),
        "nopbc": not pbc,
    }
    return info_dict
