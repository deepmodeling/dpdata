from __future__ import annotations

import numpy as np


def from_system_data(structure) -> dict:
    symbols = [site.species_string for site in structure]
    atom_names = list(structure.symbol_set)
    atom_numbs = [symbols.count(symbol) for symbol in atom_names]
    atom_types = np.array([atom_names.index(symbol) for symbol in symbols]).astype(int)
    coords = structure.cart_coords
    cells = structure.lattice.matrix

    info_dict = {
        "atom_names": atom_names,
        "atom_numbs": atom_numbs,
        "atom_types": atom_types,
        "coords": np.array([coords]),
        "cells": np.array([cells]),
    }
    return info_dict
