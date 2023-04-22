import numpy as np

from dpdata.format import Format
from dpdata.psi4.output import read_psi4_output
from dpdata.unit import EnergyConversion, ForceConversion, LengthConversion

length_convert = LengthConversion("bohr", "angstrom").value()
energy_convert = EnergyConversion("hartree", "eV").value()
force_convert = ForceConversion("hartree/bohr", "eV/angstrom").value()


@Format.register("psi4/out")
class PSI4OutFormat(Format):
    """Psi4 output.

    Note that both the energy and the gradient should be
    printed into the output file.
    """

    def from_labeled_system(self, file_name: str, **kwargs) -> dict:
        """Read from Psi4 output.

        Parameters
        ----------
        file_name : str
            file name
        **kwargs
            keyword arguments

        Returns
        -------
        dict
            system data
        """
        symbols, coord, energy, forces = read_psi4_output(file_name)

        atom_names, atom_types, atom_numbs = np.unique(
            symbols, return_inverse=True, return_counts=True
        )
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
