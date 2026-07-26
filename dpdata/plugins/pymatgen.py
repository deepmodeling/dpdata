from __future__ import annotations

import numpy as np

import dpdata.formats.pymatgen.molecule
import dpdata.formats.pymatgen.structure
from dpdata.format import Format


@Format.register("pymatgen/structure")
class PyMatgenStructureFormat(Format):
    """In-memory pymatgen ``Structure`` objects for periodic systems.

    `pymatgen <https://pymatgen.org/>`_ (Python Materials Genomics) is a
    robust, open-source Python library for materials analysis.

    This adapter converts without writing a file and requires the optional
    ``pymatgen`` dependency. Writing a multi-frame System returns one
    ``Structure`` object per frame.
    """

    def from_system(self, structure, **kwargs) -> dict:
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
        return dpdata.formats.pymatgen.structure.from_system_data(structure)

    def to_system(self, data, **kwargs):
        """Convert System frames to pymatgen Structure objects.

        Parameters
        ----------
        data : dict
            Periodic System data.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        list[pymatgen.core.Structure]
            One structure per frame.
        """
        structures = []
        try:
            from pymatgen.core import Lattice, Structure
        except ModuleNotFoundError as e:
            raise ImportError("No module pymatgen.Structure") from e

        species = [data["atom_names"][tt] for tt in data["atom_types"]]
        pbc = not (data.get("nopbc", False))
        for ii in range(data["coords"].shape[0]):
            structure = Structure(
                Lattice(data["cells"][ii], pbc=[pbc] * 3),
                species,
                data["coords"][ii],
                coords_are_cartesian=True,
            )
            structures.append(structure)
        return structures


@Format.register("pymatgen/molecule")
class PyMatgenMoleculeFormat(Format):
    """In-memory pymatgen ``Molecule`` objects for nonperiodic systems.

    `pymatgen <https://pymatgen.org/>`_ (Python Materials Genomics) is a
    robust, open-source Python library for materials analysis.

    Periodic boundary conditions are removed during conversion. The optional
    ``pymatgen`` dependency is required, and writing returns one ``Molecule``
    object per frame.
    """

    @Format.post("remove_pbc")
    def from_system(self, file_name, **kwargs):
        """Convert a pymatgen Molecule into System data.

        Parameters
        ----------
        file_name : pymatgen.core.Molecule
            In-memory molecule to convert.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        dict
            Nonperiodic system data.
        """
        try:
            from pymatgen.core import Molecule  # noqa: F401
        except ModuleNotFoundError as e:
            raise ImportError("No module pymatgen.Molecule") from e

        return dpdata.formats.pymatgen.molecule.to_system_data(file_name)

    def to_system(self, data, **kwargs):
        """Convert System frames to pymatgen Molecule objects.

        Parameters
        ----------
        data : dict
            System data. Periodic boundary conditions are removed.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        list[pymatgen.core.Molecule]
            One molecule per frame.
        """
        molecules = []
        try:
            from pymatgen.core import Molecule
        except ModuleNotFoundError as e:
            raise ImportError("No module pymatgen.Molecule") from e

        species = [data["atom_names"][tt] for tt in data["atom_types"]]
        data = dpdata.system.remove_pbc(data)
        for ii in range(np.array(data["coords"]).shape[0]):
            molecule = Molecule(species, data["coords"][ii])
            molecules.append(molecule)
        return molecules


@Format.register("pymatgen/computedstructureentry")
@Format.register_to("to_pymatgen_ComputedStructureEntry")
class PyMatgenCSEFormat(Format):
    """In-memory pymatgen ``ComputedStructureEntry`` objects.

    `pymatgen <https://pymatgen.org/>`_ (Python Materials Genomics) is a
    robust, open-source Python library for materials analysis.

    This write-only labeled adapter creates one entry per frame and places
    forces and virials in the entry's ``data`` mapping. The optional
    ``pymatgen`` dependency is required.
    """

    def to_labeled_system(self, data, *args, **kwargs):
        """Convert labeled frames to pymatgen ComputedStructureEntry objects.

        Parameters
        ----------
        data : dict
            LabeledSystem data containing energy, forces, and virials.
        *args : list
            Additional positional arguments accepted for API compatibility.
        **kwargs : dict
            Additional keyword arguments accepted for API compatibility.

        Returns
        -------
        list[pymatgen.entries.computed_entries.ComputedStructureEntry]
            One computed entry per frame.
        """
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
