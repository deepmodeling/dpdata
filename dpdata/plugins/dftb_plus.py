import numpy as np
from dpdata.format import Format
from dpdata.dftb_plus.inout import read_dftb_plus
from dpdata.unit import LengthConversion, EnergyConversion, ForceConversion


length_convert = LengthConversion("bohr", "angstrom").value()
energy_convert = EnergyConversion("hartree", "eV").value()
force_convert = ForceConversion("hartree/bohr", "eV/angstrom").value()


@Format.register("dftb_plus")
class DFTBplusFormat(Format):
    def from_labeled_system(self, file_name, **kwargs):
        symbols, coord, energy, forces = read_dftb_plus(file_name)

        atom_names, atom_types, atom_numbs = np.unique(symbols, return_inverse=True, return_counts=True)
        natoms = coord.shape[0]

        return {
            "atom_types": atom_types,
            "atom_names": list(atom_names),
            "atom_numbs": list(atom_numbs),
            "coords": (coord * length_convert).reshape((1, natoms, 3)),
            "energies": np.array([energy * energy_convert]),
            "forces": (forces * force_convert).reshape((1, natoms, 3)),
            "cells": np.zeros((1, 3, 3)),
            "orig": np.zeros(3),
            "nopbc": True,
        }
