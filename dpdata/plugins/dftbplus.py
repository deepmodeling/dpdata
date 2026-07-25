from __future__ import annotations

import numpy as np

from dpdata.format import Format
from dpdata.formats.dftbplus.output import read_dftb_plus
from dpdata.unit import EnergyConversion, ForceConversion

energy_convert = EnergyConversion("hartree", "eV").value()
force_convert = ForceConversion("hartree/bohr", "eV/angstrom").value()


@Format.register("dftbplus")
class DFTBplusFormat(Format):
    """DFTB+ input/output pair for one labeled molecular configuration.

    Pass a tuple containing the DFTB+ input geometry file and output result
    file. The reader combines symbols and coordinates from the input with the
    energy and forces from the output. See the `DFTB+ documentation
    <https://dftbplus.org/documentation>`_ for the underlying file formats.
    """

    def from_labeled_system(self, file_paths, **kwargs):
        """Reads system information from the given DFTB+ file paths.

        Parameters
        ----------
        file_paths : tuple
            A tuple containing the input and output file paths.
            - Input file (file_in): Contains information about symbols and coord.
            - Output file (file_out): Contains information about energy and force.
        **kwargs : dict
            other parameters

        """
        file_in, file_out = file_paths
        symbols, coord, energy, forces = read_dftb_plus(file_in, file_out)
        last_occurrence = {v: i for i, v in enumerate(symbols)}
        atom_names = np.array(sorted(np.unique(symbols), key=last_occurrence.get))
        atom_types = np.array([np.where(atom_names == s)[0][0] for s in symbols])
        atom_numbs = np.array([np.sum(atom_types == i) for i in range(len(atom_names))])
        natoms = coord.shape[0]

        return {
            "atom_types": atom_types,
            "atom_names": list(atom_names),
            "atom_numbs": list(atom_numbs),
            "coords": coord.reshape((1, natoms, 3)),
            "energies": np.array([energy * energy_convert]),
            "forces": (forces * force_convert).reshape((1, natoms, 3)),
            "cells": np.zeros((1, 3, 3)),
            "orig": np.zeros(3),
            "nopbc": True,
        }
