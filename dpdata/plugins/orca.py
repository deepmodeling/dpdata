from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from dpdata.format import Format
from dpdata.formats.orca.output import read_orca_sp_output
from dpdata.unit import EnergyConversion, ForceConversion

if TYPE_CHECKING:
    from dpdata.utils import FileType

energy_convert = EnergyConversion("hartree", "eV").value()
force_convert = ForceConversion("hartree/bohr", "eV/angstrom").value()


@Format.register("orca/spout")
class ORCASPOutFormat(Format):
    """ORCA single-point energy and gradient output.

    `ORCA <https://orcaforum.kofo.mpg.de/>`_ is an ab initio quantum chemistry
    program package for molecular calculations.

    The reader creates one nonperiodic :class:`dpdata.LabeledSystem` frame
    containing atomic coordinates, the total energy, and forces converted from
    the Cartesian gradient. Both energy and gradient sections must be present
    in the text output.
    """

    def from_labeled_system(self, file_name: FileType, **kwargs) -> dict:
        """Read from ORCA single point energy output.

        Parameters
        ----------
        file_name : FileType
            file name
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        dict
            system data
        """
        symbols, coord, energy, forces = read_orca_sp_output(file_name)

        atom_names, atom_types, atom_numbs = np.unique(
            symbols, return_inverse=True, return_counts=True
        )
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
