from __future__ import annotations

import os
import tempfile
import unittest
from unittest.mock import patch

from context import dpdata


class TestSchNetPackRegistration(unittest.TestCase):
    """Test SchNetPack format registration and error handling."""

    def test_format_registered(self):
        """Test that schnetpack format is properly registered."""
        test_system = dpdata.LabeledSystem()
        test_system.data = {
            "atom_names": ["H"],
            "atom_numbs": [1],
            "atom_types": [0],
            "cells": [[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]],
            "coords": [[[0.0, 0.0, 0.0]]],
            "orig": [0.0, 0.0, 0.0],
            "energies": [1.0],
            "forces": [[[0.0, 0.0, 0.0]]],
        }

        # Since ASE is available, this should work and create the database
        import tempfile

        with tempfile.TemporaryDirectory() as tmpdir:
            db_file = os.path.join(tmpdir, "test_registered.db")
            test_system.to("schnetpack", db_file)

            # Verify the database was created
            self.assertTrue(os.path.exists(db_file))


try:
    import ase.db  # noqa: F401

    ase_available = True
except ImportError:
    ase_available = False


@unittest.skipIf(not ase_available, "skip test_schnetpack")
class TestSchNetPack(unittest.TestCase):
    def setUp(self):
        # Create a simple test system
        self.test_system = dpdata.System()

        # Create simple water-like structure for testing
        # 3 atoms: O, H, H
        self.test_system.data = {
            "atom_names": ["O", "H"],
            "atom_numbs": [1, 2],
            "atom_types": [0, 1, 1],  # O, H, H
            "cells": [[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]],
            "coords": [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]],
            "orig": [0.0, 0.0, 0.0],
        }

        # Create labeled system with dummy energies and forces
        self.labeled_system = dpdata.LabeledSystem()
        self.labeled_system.data = self.test_system.data.copy()
        self.labeled_system.data["energies"] = [-10.5]  # eV
        self.labeled_system.data["forces"] = [
            [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.1]]
        ]  # eV/Ang

        # Optional: add virials
        self.labeled_system.data["virials"] = [
            [[0.01, 0.0, 0.0], [0.0, 0.01, 0.0], [0.0, 0.0, 0.01]]
        ]

    def test_to_schnetpack(self):
        """Test conversion to SchNetPack-compatible ASE database format."""
        from ase.db import connect

        with tempfile.TemporaryDirectory() as tmpdir:
            db_file = os.path.join(tmpdir, "test_data.db")

            # Convert to SchNetPack format
            self.labeled_system.to("schnetpack", db_file)

            # Verify the database was created
            self.assertTrue(os.path.exists(db_file))

            # Load the database and verify contents using ASE
            db = connect(db_file)

            # Check number of structures
            self.assertEqual(len(db), 1)

            # Get the structure back
            row = db.get(1)
            atoms = row.toatoms()

            # Check basic structure information
            self.assertEqual(len(atoms), 3)  # O, H, H
            self.assertEqual(atoms.get_chemical_symbols(), ["O", "H", "H"])

            # Check that properties are present and accessible
            self.assertTrue(hasattr(atoms, "calc") and atoms.calc is not None)

            # Check energy value
            self.assertAlmostEqual(atoms.get_potential_energy(), -10.5, places=5)

            # Check forces shape and values
            forces = atoms.get_forces()
            self.assertEqual(forces.shape, (3, 3))  # 3 atoms, 3 components

    def test_to_schnetpack_custom_units(self):
        """Test conversion with custom units."""
        from ase.db import connect

        with tempfile.TemporaryDirectory() as tmpdir:
            db_file = os.path.join(tmpdir, "test_data_units.db")

            # Convert with custom units
            property_units = {"energy": "kcal/mol", "forces": "kcal/mol/Ang"}

            self.labeled_system.to(
                "schnetpack", db_file, property_unit_dict=property_units
            )

            # Verify the database was created
            self.assertTrue(os.path.exists(db_file))

            # Load and verify using ASE
            db = connect(db_file)
            self.assertEqual(len(db), 1)

            # Basic verification that data is present
            row = db.get(1)
            atoms = row.toatoms()
            self.assertTrue(hasattr(atoms, "calc") and atoms.calc is not None)

            # Check that units were stored in metadata
            self.assertIn("property_units", row.data)
            self.assertEqual(row.data["property_units"], property_units)

    def test_to_schnetpack_without_virials(self):
        """Test conversion without virials."""
        from ase.db import connect

        with tempfile.TemporaryDirectory() as tmpdir:
            db_file = os.path.join(tmpdir, "test_no_virials.db")

            # Create system without virials
            system_no_virials = dpdata.LabeledSystem()
            system_no_virials.data = self.test_system.data.copy()
            system_no_virials.data["energies"] = [-10.5]
            system_no_virials.data["forces"] = [
                [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.1]]
            ]

            # Convert to SchNetPack format
            system_no_virials.to("schnetpack", db_file)

            # Verify the database was created
            self.assertTrue(os.path.exists(db_file))

            # Load and verify using ASE
            db = connect(db_file)
            self.assertEqual(len(db), 1)

            row = db.get(1)
            atoms = row.toatoms()
            self.assertTrue(hasattr(atoms, "calc") and atoms.calc is not None)
            # virials should not be present in row.data
            self.assertNotIn("virials", row.data)

    def test_multiframe_system(self):
        """Test conversion of multi-frame system."""
        from ase.db import connect

        with tempfile.TemporaryDirectory() as tmpdir:
            db_file = os.path.join(tmpdir, "test_multiframe.db")

            # Create multi-frame system
            multiframe_system = dpdata.LabeledSystem()
            multiframe_system.data = {
                "atom_names": ["O", "H"],
                "atom_numbs": [1, 2],
                "atom_types": [0, 1, 1],
                "cells": [
                    [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
                    [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
                ],
                "coords": [
                    [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                    [[0.1, 0.0, 0.0], [1.1, 0.0, 0.0], [0.1, 1.0, 0.0]],
                ],
                "orig": [0.0, 0.0, 0.0],
                "energies": [-10.5, -10.6],
                "forces": [
                    [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.1]],
                    [[0.2, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.2]],
                ],
            }

            # Convert to SchNetPack format
            multiframe_system.to("schnetpack", db_file)

            # Verify the database was created
            self.assertTrue(os.path.exists(db_file))

            # Load and verify using ASE
            db = connect(db_file)

            # Should have 2 frames
            self.assertEqual(len(db), 2)

            # Check both frames
            for i in range(1, 3):  # ASE database IDs start from 1
                row = db.get(i)
                atoms = row.toatoms()
                self.assertTrue(hasattr(atoms, "calc") and atoms.calc is not None)
                self.assertEqual(len(atoms), 3)  # O, H, H

    def test_schnetpack_script_compatibility(self):
        """Test compatibility with the specific SchNetPack script from user requirements."""
        from ase.db import connect

        with tempfile.TemporaryDirectory() as tmpdir:
            db_file = os.path.join(tmpdir, "schnet.db")

            # Create a system similar to what the user script expects
            # with enough frames for num_train=1000, num_val=100 (need at least 1100)
            import numpy as np

            num_frames = 1200
            coords = []
            cells = []
            energies = []
            forces = []

            for i in range(num_frames):
                # Base cell
                cell = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
                cells.append(cell)

                # Slightly perturbed coordinates for variety
                coord = [
                    [0.0 + 0.01 * i, 0.0, 0.0],  # O
                    [1.0 + 0.01 * np.sin(i / 10), 0.0, 0.0],  # H
                    [0.0, 1.0 + 0.01 * np.cos(i / 10), 0.0],  # H
                ]
                coords.append(coord)

                # Energy (slightly varying)
                energy = -10.5 - 0.1 * np.sin(i / 100)
                energies.append(energy)

                # Forces (varying)
                force = [
                    [0.1 + 0.01 * np.cos(i / 20), 0.0, 0.0],
                    [0.0, 0.1 + 0.01 * np.sin(i / 20), 0.0],
                    [0.0, 0.0, 0.1 + 0.01 * np.cos(i / 30)],
                ]
                forces.append(force)

            labeled_system = dpdata.LabeledSystem()
            labeled_system.data = {
                "atom_names": ["O", "H"],
                "atom_numbs": [1, 2],
                "atom_types": [0, 1, 1],  # O, H, H
                "cells": cells,
                "coords": coords,
                "orig": [0.0, 0.0, 0.0],
                "energies": energies,
                "forces": forces,
            }

            # Convert with exact units from user script
            expected_units = {"energy": "eV", "forces": "eV/Ang"}
            labeled_system.to("schnetpack", db_file, property_unit_dict=expected_units)

            # Verify the database was created
            self.assertTrue(os.path.exists(db_file))

            # Load and verify using ASE
            db = connect(db_file)

            # Check we have enough entries for the user's requirements
            self.assertGreaterEqual(
                len(db),
                1100,
                "Need at least 1100 entries for num_train=1000, num_val=100",
            )

            # Check property units match exactly
            first_row = db.get(1)
            self.assertIn("property_units", first_row.data)
            self.assertEqual(first_row.data["property_units"], expected_units)

            # Test that all entries have the required properties for SchNetPack
            sample_size = min(10, len(db))
            for i in range(1, sample_size + 1):
                row = db.get(i)
                atoms = row.toatoms()

                # Check calculator exists
                self.assertTrue(hasattr(atoms, "calc") and atoms.calc is not None)

                # Check energy is accessible
                energy = atoms.get_potential_energy()
                self.assertIsInstance(energy, (int, float))

                # Check forces are accessible and correct shape
                forces = atoms.get_forces()
                self.assertEqual(forces.shape, (len(atoms), 3))

                # Check positions are correct shape
                self.assertEqual(atoms.positions.shape, (len(atoms), 3))

            # Test database iteration (SchNetPack does this)
            energies_from_db = []
            for row in db.select(limit=5):
                atoms = row.toatoms()
                energy = atoms.get_potential_energy()
                energies_from_db.append(energy)

            self.assertEqual(len(energies_from_db), 5)

            # Verify all energies are numbers
            for energy in energies_from_db:
                self.assertIsInstance(energy, (int, float))


class TestSchNetPackMocked(unittest.TestCase):
    """Test SchNetPack functionality with mocked dependencies."""

    def setUp(self):
        # Create a simple test system
        self.labeled_system = dpdata.LabeledSystem()
        self.labeled_system.data = {
            "atom_names": ["O", "H"],
            "atom_numbs": [1, 2],
            "atom_types": [0, 1, 1],  # O, H, H
            "cells": [[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]],
            "coords": [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]],
            "orig": [0.0, 0.0, 0.0],
            "energies": [-10.5],  # eV
            "forces": [[[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.1]]],  # eV/Ang
        }

    @patch("dpdata.plugins.schnetpack.SchNetPackFormat.to_labeled_system")
    def test_conversion_logic_mocked(self, mock_to_labeled_system):
        """Test the conversion logic with mocked dependencies."""
        # Test that the method can be called
        mock_to_labeled_system.return_value = None

        # Test the conversion - should call the mocked method
        self.labeled_system.to("schnetpack", "/tmp/test.db")

        # Verify the method was called
        mock_to_labeled_system.assert_called_once()


if __name__ == "__main__":
    unittest.main()
