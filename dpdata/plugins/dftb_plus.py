import numpy as np

from dpdata.dftb_plus.inout import read_dftb_plus
from dpdata.format import Format
from dpdata.unit import EnergyConversion, ForceConversion, LengthConversion

length_convert = LengthConversion("bohr", "angstrom").value()
energy_convert = EnergyConversion("hartree", "eV").value()
force_convert = ForceConversion("hartree/bohr", "eV/angstrom").value()


@Format.register("dftb_plus")
class DFTBplusFormat(Format):
    def from_labeled_system(self, file_name, **kwargs):
        symbols, coord, energy, forces = read_dftb_plus(file_name)
        last_occurrence = {v: i for i, v in enumerate(symbols)}
        atom_names = np.array(sorted(np.unique(symbols), key=last_occurrence.get))
        atom_types = np.array([np.where(atom_names == s)[0][0] for s in symbols])
        atom_numbs = np.array([np.sum(atom_types == i) for i in range(len(atom_names))])
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
