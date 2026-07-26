from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from dpdata.format import Format
from dpdata.formats.psi4.input import write_psi4_input
from dpdata.formats.psi4.output import read_psi4_output
from dpdata.unit import EnergyConversion, ForceConversion
from dpdata.utils import open_file

if TYPE_CHECKING:
    from dpdata.utils import FileType

energy_convert = EnergyConversion("hartree", "eV").value()
force_convert = ForceConversion("hartree/bohr", "eV/angstrom").value()


@Format.register("psi4/out")
class PSI4OutFormat(Format):
    """Psi4 energy and gradient output.

    `Psi4 <https://psicode.org/>`_ is an open-source quantum chemistry program
    package.

    The reader creates one nonperiodic :class:`dpdata.LabeledSystem` frame
    containing atomic coordinates, the total energy, and forces converted from
    the Cartesian gradient. Both energy and gradient sections must be present
    in the text output.
    """

    def from_labeled_system(self, file_name: FileType, **kwargs) -> dict:
        """Read from Psi4 output.

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
        symbols, coord, energy, forces = read_psi4_output(file_name)

        atom_names, atom_types, atom_numbs = np.unique(
            symbols, return_inverse=True, return_counts=True
        )
        natoms = coord.shape[0]

        return {
            "atom_types": atom_types,
            "atom_names": list(atom_names),
            "atom_numbs": list(atom_numbs),
            "coords": (coord).reshape((1, natoms, 3)),
            "energies": np.array([energy * energy_convert]),
            "forces": (forces * force_convert).reshape((1, natoms, 3)),
            "cells": np.zeros((1, 3, 3)),
            "orig": np.zeros(3),
            "nopbc": True,
        }


@Format.register("psi4/inp")
class PSI4InputFormat(Format):
    """Psi4 input file for a single molecular configuration.

    `Psi4 <https://psicode.org/>`_ is an open-source quantum chemistry program.

    The writer emits the molecule, charge, multiplicity, requested electronic-
    structure method, and basis set. Psi4 input is write-only in dpdata; use
    ``psi4/out`` to load calculated energies and gradients.
    """

    def to_system(
        self,
        data: dict,
        file_name: FileType,
        method: str,
        basis: str,
        charge: int = 0,
        multiplicity: int = 1,
        frame_idx=0,
        **kwargs,
    ):
        """Write PSI4 input.

        Parameters
        ----------
        data : dict
            system data
        file_name : str
            file name
        method : str
            computational method
        basis : str
            basis set; see https://psicode.org/psi4manual/master/basissets_tables.html
        charge : int, default=0
            charge of system
        multiplicity : int, default=1
            multiplicity of system
        frame_idx : int, default=0
            The index of the frame to dump
        **kwargs : dict
            Additional format arguments accepted for API compatibility.
        """
        types = np.array(data["atom_names"])[data["atom_types"]]
        with open_file(file_name, "w") as fout:
            fout.write(
                write_psi4_input(
                    types=types,
                    coords=data["coords"][frame_idx],
                    method=method,
                    basis=basis,
                    charge=charge,
                    multiplicity=multiplicity,
                )
            )
