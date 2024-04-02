import numpy as np

import dpdata.pymatgen.molecule
from dpdata.format import Format
from pymatgen.core import Structure

@Format.register("pymatgen/structure")
class PyMatgenStructureFormat(Format):

    def from_system(self, structure: Structure, **kwargs) -> dict:
        """Convert pymatgen.core.Structure to System.

        Parameters
        ----------
        structure : pymatgen.core.Structure
            a Pymatgen Structure, containing a structure
        **kwargs : dict
            other parameters

        Returns
        -------
        dict
            data dict
        """
        symbols = [site.species_string for site in structure]
        atom_names = list(structure.symbol_set)
        atom_numbs = [symbols.count(symbol) for symbol in atom_names]
        atom_types = np.array([atom_names.index(symbol) for symbol in symbols]).astype(int)
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


    def to_system(self, data, **kwargs):
        """Convert System to Pymatgen Structure obj."""
        structures = []
        try:
            from pymatgen.core import Structure
        except ModuleNotFoundError as e:
            raise ImportError("No module pymatgen.Structure") from e

        species = []
        for name, numb in zip(data["atom_names"], data["atom_numbs"]):
            species.extend([name] * numb)
        for ii in range(data["coords"].shape[0]):
            structure = Structure(
                data["cells"][ii],
                species,
                data["coords"][ii],
                coords_are_cartesian=True,
            )
            structures.append(structure)
        return structures


@Format.register("pymatgen/molecule")
class PyMatgenMoleculeFormat(Format):
    @Format.post("remove_pbc")
    def from_system(self, file_name, **kwargs):
        try:
            from pymatgen.core import Molecule  # noqa: F401
        except ModuleNotFoundError as e:
            raise ImportError("No module pymatgen.Molecule") from e

        return dpdata.pymatgen.molecule.to_system_data(file_name)

    def to_system(self, data, **kwargs):
        """Convert System to Pymatgen Molecule obj."""
        molecules = []
        try:
            from pymatgen.core import Molecule
        except ModuleNotFoundError as e:
            raise ImportError("No module pymatgen.Molecule") from e

        species = []
        for name, numb in zip(data["atom_names"], data["atom_numbs"]):
            species.extend([name] * numb)
        data = dpdata.system.remove_pbc(data)
        for ii in range(np.array(data["coords"]).shape[0]):
            molecule = Molecule(species, data["coords"][ii])
            molecules.append(molecule)
        return molecules


@Format.register("pymatgen/computedstructureentry")
@Format.register_to("to_pymatgen_ComputedStructureEntry")
class PyMatgenCSEFormat(Format):
    def to_labeled_system(self, data, *args, **kwargs):
        """Convert System to Pymagen ComputedStructureEntry obj."""
        try:
            from pymatgen.entries.computed_entries import ComputedStructureEntry
        except ModuleNotFoundError as e:
            raise ImportError(
                "No module ComputedStructureEntry in pymatgen.entries.computed_entries"
            ) from e

        entries = []

        for ii, structure in enumerate(PyMatgenStructureFormat().to_system(data)):
            energy = data["energies"][ii]
            csedata = {"forces": data["forces"][ii], "virials": data["virials"][ii]}

            entry = ComputedStructureEntry(structure, energy, data=csedata)
            entries.append(entry)
        return entries
