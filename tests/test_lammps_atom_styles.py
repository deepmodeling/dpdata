from __future__ import annotations

import os
import unittest

import numpy as np
from context import dpdata


class TestLammpsAtomStyles(unittest.TestCase):
    """Test support for different LAMMPS atom styles."""

    def setUp(self):
        """Set up test fixtures."""
        # Create test data files for different atom styles
        self.test_files = {}

        # Full style test file
        full_content = """# LAMMPS data file - full style
2 atoms
2 atom types
0.0 2.5243712 xlo xhi
0.0 2.0430257 ylo yhi
0.0 2.2254033 zlo zhi
1.2621856 1.2874292 0.7485898 xy xz yz

Atoms # full

1 1 1 -0.8476 0.0 0.0 0.0
2 1 2 0.4238 1.2621856 0.7018028 0.5513885"""
        self.test_files["full"] = "/tmp/test_full_style.lmp"
        with open(self.test_files["full"], "w") as f:
            f.write(full_content)

        # Charge style test file
        charge_content = """# LAMMPS data file - charge style
2 atoms
2 atom types
0.0 2.5243712 xlo xhi
0.0 2.0430257 ylo yhi
0.0 2.2254033 zlo zhi
1.2621856 1.2874292 0.7485898 xy xz yz

Atoms # charge

1 1 -0.8476 0.0 0.0 0.0
2 2 0.4238 1.2621856 0.7018028 0.5513885"""
        self.test_files["charge"] = "/tmp/test_charge_style.lmp"
        with open(self.test_files["charge"], "w") as f:
            f.write(charge_content)

        # Bond style test file
        bond_content = """# LAMMPS data file - bond style
2 atoms
2 atom types
0.0 2.5243712 xlo xhi
0.0 2.0430257 ylo yhi
0.0 2.2254033 zlo zhi
1.2621856 1.2874292 0.7485898 xy xz yz

Atoms # bond

1 1 1 0.0 0.0 0.0
2 1 2 1.2621856 0.7018028 0.5513885"""
        self.test_files["bond"] = "/tmp/test_bond_style.lmp"
        with open(self.test_files["bond"], "w") as f:
            f.write(bond_content)

    def tearDown(self):
        """Clean up test files."""
        for file_path in self.test_files.values():
            if os.path.exists(file_path):
                os.remove(file_path)

    def test_atomic_style_backward_compatibility(self):
        """Test that atomic style still works (backward compatibility)."""
        system = dpdata.System(os.path.join("poscars", "conf.lmp"), type_map=["O", "H"])
        self.assertEqual(len(system["atom_types"]), 2)
        self.assertEqual(system["atom_types"][0], 0)  # O
        self.assertEqual(system["atom_types"][1], 1)  # H

    def test_full_style_parsing(self):
        """Test parsing of full style LAMMPS data file."""
        system = dpdata.System(
            self.test_files["full"],
            fmt="lammps/lmp",
            type_map=["O", "H"],
            atom_style="full",
        )

        # Check basic structure
        self.assertEqual(len(system["atom_types"]), 2)
        self.assertEqual(system["atom_types"][0], 0)  # type 1 -> O
        self.assertEqual(system["atom_types"][1], 1)  # type 2 -> H

        # Check coordinates
        expected_coords = np.array([[0.0, 0.0, 0.0], [1.2621856, 0.7018028, 0.5513885]])
        np.testing.assert_allclose(system["coords"][0], expected_coords, atol=1e-6)

        # Check charges are present
        self.assertIn("charges", system.data)
        expected_charges = np.array([-0.8476, 0.4238])
        np.testing.assert_allclose(system["charges"][0], expected_charges, atol=1e-6)

    def test_charge_style_parsing(self):
        """Test parsing of charge style LAMMPS data file."""
        system = dpdata.System(
            self.test_files["charge"],
            fmt="lammps/lmp",
            type_map=["O", "H"],
            atom_style="charge",
        )

        # Check basic structure
        self.assertEqual(len(system["atom_types"]), 2)
        self.assertEqual(system["atom_types"][0], 0)  # type 1 -> O
        self.assertEqual(system["atom_types"][1], 1)  # type 2 -> H

        # Check coordinates
        expected_coords = np.array([[0.0, 0.0, 0.0], [1.2621856, 0.7018028, 0.5513885]])
        np.testing.assert_allclose(system["coords"][0], expected_coords, atol=1e-6)

        # Check charges are present
        self.assertIn("charges", system.data)
        expected_charges = np.array([-0.8476, 0.4238])
        np.testing.assert_allclose(system["charges"][0], expected_charges, atol=1e-6)

    def test_bond_style_parsing(self):
        """Test parsing of bond style LAMMPS data file."""
        system = dpdata.System(
            self.test_files["bond"],
            fmt="lammps/lmp",
            type_map=["O", "H"],
            atom_style="bond",
        )

        # Check basic structure
        self.assertEqual(len(system["atom_types"]), 2)
        self.assertEqual(system["atom_types"][0], 0)  # type 1 -> O
        self.assertEqual(system["atom_types"][1], 1)  # type 2 -> H

        # Check coordinates
        expected_coords = np.array([[0.0, 0.0, 0.0], [1.2621856, 0.7018028, 0.5513885]])
        np.testing.assert_allclose(system["coords"][0], expected_coords, atol=1e-6)

        # Bond style should not have charges
        self.assertNotIn("charges", system.data)

    def test_unsupported_atom_style(self):
        """Test that unsupported atom styles raise appropriate errors."""
        with self.assertRaises(ValueError) as context:
            dpdata.System(
                self.test_files["bond"],
                fmt="lammps/lmp",
                type_map=["O", "H"],
                atom_style="unsupported_style",
            )

        self.assertIn("Unsupported atom style", str(context.exception))

    def test_default_atomic_style(self):
        """Test that default behavior is atomic style."""
        # Test using existing atomic style file
        system1 = dpdata.System(
            os.path.join("poscars", "conf.lmp"), type_map=["O", "H"]
        )
        system2 = dpdata.System(
            os.path.join("poscars", "conf.lmp"),
            type_map=["O", "H"],
            atom_style="atomic",
        )

        # Should be identical
        np.testing.assert_array_equal(system1["coords"], system2["coords"])
        np.testing.assert_array_equal(system1["atom_types"], system2["atom_types"])


if __name__ == "__main__":
    unittest.main()
