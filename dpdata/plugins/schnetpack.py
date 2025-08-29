from __future__ import annotations

from dpdata.format import Format


@Format.register("schnetpack")
class SchNetPackFormat(Format):
    """Format for SchNetPack ASE database.

    SchNetPack uses ASE database format internally for storing atomic structures
    and their properties. This format converts dpdata LabeledSystem to SchNetPack's
    ASE database format.

    For more information, see:
    https://schnetpack.readthedocs.io/en/latest/tutorials/tutorial_01_preparing_data.html
    """

    def to_labeled_system(
        self,
        data: dict,
        file_name: str = "schnetpack_data.db",
        distance_unit: str = "Ang",
        property_unit_dict: dict | None = None,
        **kwargs,
    ) -> None:
        """Convert dpdata LabeledSystem to SchNetPack ASE database format.

        Parameters
        ----------
        data : dict
            dpdata LabeledSystem data dictionary
        file_name : str, optional
            Path to the output SchNetPack database file, by default "schnetpack_data.db"
        distance_unit : str, optional
            Unit for distances, by default "Ang"
        property_unit_dict : dict, optional
            Dictionary mapping property names to their units.
            If None, defaults to {"energy": "eV", "forces": "eV/Ang"}
        **kwargs : dict
            Additional keyword arguments

        Raises
        ------
        ImportError
            If ASE or SchNetPack are not available
        """
        try:
            from ase import Atoms
            from schnetpack.data import ASEAtomsData
        except ImportError as e:
            raise ImportError(
                "ASE and SchNetPack are required for schnetpack format. "
                "Install with: pip install ase schnetpack"
            ) from e

        # Set default units if not provided
        if property_unit_dict is None:
            property_unit_dict = {"energy": "eV", "forces": "eV/Ang"}

        # Convert dpdata to list of ASE Atoms and property list
        atoms_list = []
        property_list = []

        species = [data["atom_names"][tt] for tt in data["atom_types"]]
        nframes = data["coords"].shape[0]

        for frame_idx in range(nframes):
            # Create ASE Atoms object for this frame
            atoms = Atoms(
                symbols=species,
                positions=data["coords"][frame_idx],
                pbc=not data.get("nopbc", False),
                cell=data["cells"][frame_idx],
            )
            atoms_list.append(atoms)

            # Create property dictionary for this frame
            properties = {}

            # Add energy
            if "energies" in data:
                properties["energy"] = float(data["energies"][frame_idx])

            # Add forces
            if "forces" in data:
                properties["forces"] = data["forces"][frame_idx]

            # Add virials if present (SchNetPack doesn't have built-in support,
            # but can be stored as additional property)
            if "virials" in data:
                properties["virials"] = data["virials"][frame_idx]

            property_list.append(properties)

        # Create SchNetPack ASE database
        dataset = ASEAtomsData.create(
            file_name,
            distance_unit=distance_unit,
            property_unit_dict=property_unit_dict,
        )

        # Add all systems to the database
        dataset.add_systems(property_list, atoms_list)

        return None
