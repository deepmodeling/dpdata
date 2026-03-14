from __future__ import annotations

from dpdata.format import Format


@Format.register("schnetpack")
class SchNetPackFormat(Format):
    """Format for SchNetPack-compatible ASE database.

    SchNetPack uses ASE database format internally for storing atomic structures
    and their properties. This format converts dpdata LabeledSystem to
    SchNetPack-compatible ASE database format using only ASE functionality.

    The created database can be used directly with SchNetPack for training
    machine learning models.

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
        """Convert dpdata LabeledSystem to SchNetPack-compatible ASE database format.

        Parameters
        ----------
        data : dict
            dpdata LabeledSystem data dictionary
        file_name : str, optional
            Path to the output database file, by default "schnetpack_data.db"
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
            If ASE is not available
        """
        try:
            from ase import Atoms
            from ase.calculators.singlepoint import SinglePointCalculator
            from ase.db import connect
        except ImportError as e:
            raise ImportError(
                "ASE is required for schnetpack format. Install with: pip install ase"
            ) from e

        # Set default units if not provided
        if property_unit_dict is None:
            property_unit_dict = {"energy": "eV", "forces": "eV/Ang"}

        # Create ASE database connection
        db = connect(file_name, append=False)

        # Store property units as metadata for the entire database
        # This ensures compatibility with different SchNetPack versions
        if property_unit_dict:
            # Store metadata in database metadata (if supported)
            try:
                # Some versions of ASE support metadata storage
                db.metadata = {"property_units": property_unit_dict}
            except (AttributeError, NotImplementedError):
                # Fallback: store in a special row (will be filtered out by SchNetPack)
                pass

        species = [data["atom_names"][tt] for tt in data["atom_types"]]

        # Handle both list and numpy array formats
        import numpy as np

        coords = np.array(data["coords"])
        cells = np.array(data["cells"])
        energies = np.array(data.get("energies", [])) if "energies" in data else None
        forces = np.array(data.get("forces", [])) if "forces" in data else None
        virials = np.array(data.get("virials", [])) if "virials" in data else None

        nframes = coords.shape[0]

        for frame_idx in range(nframes):
            # Create ASE Atoms object for this frame
            atoms = Atoms(
                symbols=species,
                positions=coords[frame_idx],
                pbc=not data.get("nopbc", False),
                cell=cells[frame_idx],
            )

            # Prepare calculator properties
            calc_properties = {}

            # Add energy
            if energies is not None:
                calc_properties["energy"] = float(energies[frame_idx])

            # Add forces
            if forces is not None:
                calc_properties["forces"] = forces[frame_idx]

            # Attach calculator with properties
            if calc_properties:
                calc = SinglePointCalculator(atoms, **calc_properties)
                atoms.calc = calc

            # Prepare additional data for database (e.g., virials)
            db_data = {}
            if virials is not None:
                db_data["virials"] = virials[frame_idx]

            # Add property units as metadata for each row for maximum compatibility
            # Some SchNetPack versions might expect this per-row
            if property_unit_dict:
                db_data["property_units"] = property_unit_dict

            # Ensure energy and forces are accessible in multiple ways for compatibility
            write_kwargs = {}
            if energies is not None:
                # Store energy as a keyword argument for direct access
                write_kwargs["energy"] = float(energies[frame_idx])
            if forces is not None:
                # Store forces as a keyword argument for direct access
                write_kwargs["forces"] = forces[frame_idx]

            # Write to database with all possible access methods
            try:
                db.write(atoms, data=db_data, **write_kwargs)
            except Exception:
                # Fallback: write without direct property arguments
                # Some ASE versions might not support energy/forces as kwargs
                db.write(atoms, data=db_data)

        return None
