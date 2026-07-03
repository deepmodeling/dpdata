from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from dpdata.format import Format
from dpdata.utils import open_file

if TYPE_CHECKING:
    from dpdata.utils import FileType

from ..unit import EnergyConversion, ForceConversion, LengthConversion

length_convert = LengthConversion("bohr", "angstrom").value()
energy_convert = EnergyConversion("hartree", "eV").value()
force_convert = ForceConversion("hartree/bohr", "eV/angstrom").value()


def match_indices(atype1, atype2):
    # Ensure atype2 is a numpy array for efficient operations
    atype2 = np.array(atype2)
    # Placeholder for matched indices
    matched_indices = []
    # Track used indices to handle duplicates
    used_indices = set()

    # Iterate over each element in atype1
    for element in atype1:
        # Find all indices of the current element in atype2
        # np.where returns a tuple, so [0] is used to access the array of indices
        indices = np.where(atype2 == element)[0]

        # Find the first unused index
        for index in indices:
            if index not in used_indices:
                # Add the index to the results and mark it as used
                matched_indices.append(index)
                used_indices.add(index)
                break  # Move to the next element in atype1

    return matched_indices


@Format.register("n2p2")
class N2P2Format(Format):
    """n2p2.

    This class support the conversion from and to the training data of n2p2 format.
    For more information about the n2p2 format, please refer to https://compphysvienna.github.io/n2p2/topics/cfg_file.html
    """

    def from_labeled_system(self, file_name: FileType, **kwargs):
        """Read from n2p2 format.

        Parameters
        ----------
        file_name : str
            file name, i.e. the first argument
        **kwargs : dict
            keyword arguments that will be passed from the method

        Returns
        -------
        data : dict
            system data, whose keys are defined in LabeledSystem.DTYPES
        """
        cells = []
        coords = []
        atypes = []
        forces = []
        energies = []
        natom0 = None
        natoms0 = None
        atom_types0 = None
        with open_file(file_name) as file:
            for line in file:
                line = line.strip()  # Remove leading/trailing whitespace
                if line.lower() == "begin":
                    current_section = []  # Start a new section
                    cell = []
                    coord = []
                    atype = []
                    force = []
                    energy = None
                elif line.lower() == "end":
                    # If we are at the end of a section, process the section
                    assert len(coord) == len(atype) == len(force), (
                        "Number of atoms, atom types, and forces must match."
                    )

                    # Check if the number of atoms is consistent across all frames
                    natom = len(coord)
                    if natom0 is None:
                        natom0 = natom
                    else:
                        assert natom == natom0, (
                            "The number of atoms in all frames must be the same."
                        )

                    # Check if the number of atoms of each type is consistent across all frames
                    atype = np.array(atype)
                    unique_dict = {element: None for element in atype}
                    unique_atypes = np.array(list(unique_dict.keys()))
                    unique_atypes_list = list(unique_atypes)
                    ntypes = len(unique_atypes)
                    natoms = [len(atype[atype == at]) for at in unique_atypes]
                    if natoms0 is None:
                        natoms0 = natoms
                    else:
                        assert natoms == natoms0, (
                            "The number of atoms of each type in all frames must be the same."
                        )
                    if atom_types0 is None:
                        atom_types0 = atype
                    atom_order = match_indices(atom_types0, atype)

                    cell = np.array(cell, dtype=float)
                    coord = np.array(coord, dtype=float)[atom_order]
                    force = np.array(force, dtype=float)[atom_order]

                    cells.append(cell)
                    coords.append(coord)
                    forces.append(force)
                    energies.append(float(energy))

                    current_section = None  # Reset for the next section
                elif current_section is not None:
                    # If we are inside a section, append the line to the current section
                    # current_section.append(line)
                    line_contents = line.split()
                    if line_contents[0] == "lattice":
                        cell.append(line_contents[1:])
                    elif line_contents[0] == "atom":
                        coord.append(line_contents[1:4])
                        atype.append(line_contents[4])
                        force.append(line_contents[7:10])
                    elif line_contents[0] == "energy":
                        energy = line_contents[1]

        atom_names = unique_atypes_list
        atom_numbs = natoms
        atom_types = np.zeros(len(atom_types0), dtype=int)
        for i in range(ntypes):
            atom_types[atom_types0 == unique_atypes_list[i]] = i

        cells = np.array(cells) * length_convert
        coords = np.array(coords) * length_convert
        forces = np.array(forces) * force_convert
        energies = np.array(energies) * energy_convert

        return {
            "atom_names": list(atom_names),
            "atom_numbs": list(atom_numbs),
            "atom_types": atom_types,
            "coords": coords,
            "cells": cells,
            "nopbc": False,
            "orig": np.zeros(3),
            "energies": energies,
            "forces": forces,
        }

    def to_labeled_system(self, data, file_name: FileType, **kwargs):
        """Write n2p2 format.

        By default, LabeledSystem.to will fallback to System.to.

        Parameters
        ----------
        data : dict
            system data, whose keys are defined in LabeledSystem.DTYPES
        file_name : str
            file name, where the data will be written
        *args : list
            arguments that will be passed from the method
        **kwargs : dict
            keyword arguments that will be passed from the method
        """
        buff = []
        nframe = len(data["energies"])
        natom = len(data["atom_types"])
        atom_names = data["atom_names"]
        for frame in range(nframe):
            coord = data["coords"][frame] / length_convert
            force = data["forces"][frame] / force_convert
            energy = data["energies"][frame] / energy_convert
            cell = data["cells"][frame] / length_convert
            atype = data["atom_types"]
            buff.append("begin")
            for i in range(3):
                buff.append(
                    f"lattice {cell[i][0]:15.6f}  {cell[i][1]:15.6f}  {cell[i][2]:15.6f}"
                )
            for i in range(natom):
                buff.append(
                    f"atom {coord[i][0]:15.6f} {coord[i][1]:15.6f} {coord[i][2]:15.6f} {atom_names[atype[i]]:>7} {0:15.6f} {0:15.6f} {force[i][0]:15.6e} {force[i][1]:15.6e} {force[i][2]:15.6e}"
                )
            buff.append(f"energy {energy:15.6f}")
            buff.append(f"charge {0:15.6f}")
            buff.append("end")
        with open_file(file_name, "w") as fp:
            fp.write("\n".join(buff))
