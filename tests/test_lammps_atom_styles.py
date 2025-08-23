from __future__ import annotations

import os
import unittest

import numpy as np
from context import dpdata


class TestLammpsAtomStyles(unittest.TestCase):
    """Test support for different LAMMPS atom styles."""

    def setUp(self) -> None:
        """Set up test fixtures."""
        # Create test data files for different atom styles
        self.test_files = {}

        # Test data configurations
        self.test_configs = {
            "full": {
                "content": """# LAMMPS data file - full style
2 atoms
2 atom types
0.0 2.5243712 xlo xhi
0.0 2.0430257 ylo yhi
0.0 2.2254033 zlo zhi
1.2621856 1.2874292 0.7485898 xy xz yz

Atoms # full

1 1 1 -0.8476 0.0 0.0 0.0
2 1 2 0.4238 1.2621856 0.7018028 0.5513885""",
                "has_charges": True,
                "expected_charges": [-0.8476, 0.4238],
                "expected_coords": [[0.0, 0.0, 0.0], [1.2621856, 0.7018028, 0.5513885]],
            },
            "charge": {
                "content": """# LAMMPS data file - charge style
2 atoms
2 atom types
0.0 2.5243712 xlo xhi
0.0 2.0430257 ylo yhi
0.0 2.2254033 zlo zhi
1.2621856 1.2874292 0.7485898 xy xz yz

Atoms # charge

1 1 -0.8476 0.0 0.0 0.0
2 2 0.4238 1.2621856 0.7018028 0.5513885""",
                "has_charges": True,
                "expected_charges": [-0.8476, 0.4238],
                "expected_coords": [[0.0, 0.0, 0.0], [1.2621856, 0.7018028, 0.5513885]],
            },
            "bond": {
                "content": """# LAMMPS data file - bond style
2 atoms
2 atom types
0.0 2.5243712 xlo xhi
0.0 2.0430257 ylo yhi
0.0 2.2254033 zlo zhi
1.2621856 1.2874292 0.7485898 xy xz yz

Atoms # bond

1 1 1 0.0 0.0 0.0
2 1 2 1.2621856 0.7018028 0.5513885""",
                "has_charges": False,
                "expected_charges": None,
                "expected_coords": [[0.0, 0.0, 0.0], [1.2621856, 0.7018028, 0.5513885]],
            },
            # Test files without style comments for heuristic detection
            "full_no_comment": {
                "content": """# LAMMPS data file - full style without comment
2 atoms
2 atom types
0.0 2.5243712 xlo xhi
0.0 2.0430257 ylo yhi
0.0 2.2254033 zlo zhi
1.2621856 1.2874292 0.7485898 xy xz yz

Atoms

1 1 1 -0.8476 0.0 0.0 0.0
2 1 2 0.4238 1.2621856 0.7018028 0.5513885""",
                "has_charges": True,
                "expected_charges": [-0.8476, 0.4238],
                "expected_coords": [[0.0, 0.0, 0.0], [1.2621856, 0.7018028, 0.5513885]],
            },
            "charge_no_comment": {
                "content": """# LAMMPS data file - charge style without comment
2 atoms
2 atom types
0.0 2.5243712 xlo xhi
0.0 2.0430257 ylo yhi
0.0 2.2254033 zlo zhi
1.2621856 1.2874292 0.7485898 xy xz yz

Atoms

1 1 -0.8476 0.0 0.0 0.0
2 2 0.4238 1.2621856 0.7018028 0.5513885""",
                "has_charges": True,
                "expected_charges": [-0.8476, 0.4238],
                "expected_coords": [[0.0, 0.0, 0.0], [1.2621856, 0.7018028, 0.5513885]],
            },
            "bond_no_comment": {
                "content": """# LAMMPS data file - bond style without comment
2 atoms
2 atom types
0.0 2.5243712 xlo xhi
0.0 2.0430257 ylo yhi
0.0 2.2254033 zlo zhi
1.2621856 1.2874292 0.7485898 xy xz yz

Atoms

1 1 1 0.0 0.0 0.0
2 1 2 1.2621856 0.7018028 0.5513885""",
                "has_charges": False,
                "expected_charges": None,
                "expected_coords": [[0.0, 0.0, 0.0], [1.2621856, 0.7018028, 0.5513885]],
            },
        }

        # Create test files
        for style, config in self.test_configs.items():
            filepath = f"/tmp/test_{style}_style.lmp"
            self.test_files[style] = filepath
            with open(filepath, "w") as f:
                f.write(config["content"])

    def tearDown(self) -> None:
        """Clean up test files."""
        for file_path in self.test_files.values():
            if os.path.exists(file_path):
                os.remove(file_path)

    def _load_system(
        self, style: str, explicit_style: str | None = None
    ) -> dpdata.System:
        """Helper method to load a system with the given style.

        Parameters
        ----------
        style : str
            The style configuration key from self.test_configs
        explicit_style : str, optional
            Explicit atom_style parameter to pass to System()

        Returns
        -------
        dpdata.System
            Loaded system
        """
        kwargs = {
            "file_name": self.test_files[style],
            "fmt": "lammps/lmp",
            "type_map": ["O", "H"],
        }
        if explicit_style is not None:
            kwargs["atom_style"] = explicit_style

        return dpdata.System(**kwargs)

    def _assert_basic_structure(self, system: dpdata.System) -> None:
        """Helper method to check basic system structure."""
        self.assertEqual(len(system["atom_types"]), 2)
        self.assertEqual(system["atom_types"][0], 0)  # type 1 -> O
        self.assertEqual(system["atom_types"][1], 1)  # type 2 -> H

    def _assert_coordinates(
        self, system: dpdata.System, expected_coords: list[list[float]]
    ) -> None:
        """Helper method to check coordinates."""
        np.testing.assert_allclose(system["coords"][0], expected_coords, atol=1e-6)

    def _assert_charges(
        self, system: dpdata.System, expected_charges: list[float] | None
    ) -> None:
        """Helper method to check charges."""
        if expected_charges is not None:
            self.assertIn("charges", system.data)
            np.testing.assert_allclose(
                system["charges"][0], expected_charges, atol=1e-6
            )
        else:
            self.assertNotIn("charges", system.data)

    def _test_style_parsing(
        self, style_key: str, explicit_style: str | None = None
    ) -> None:
        """Generic helper method to test style parsing.

        Parameters
        ----------
        style_key : str
            Key from self.test_configs to test
        explicit_style : str, optional
            Explicit atom_style to pass (for testing backward compatibility)
        """
        config = self.test_configs[style_key]
        system = self._load_system(style_key, explicit_style)

        # Check basic structure
        self._assert_basic_structure(system)

        # Check coordinates
        self._assert_coordinates(system, config["expected_coords"])

        # Check charges
        self._assert_charges(system, config["expected_charges"])

    def test_atomic_style_backward_compatibility(self) -> None:
        """Test that atomic style still works (backward compatibility)."""
        system = dpdata.System(os.path.join("poscars", "conf.lmp"), type_map=["O", "H"])
        self.assertEqual(len(system["atom_types"]), 2)
        self.assertEqual(system["atom_types"][0], 0)  # O
        self.assertEqual(system["atom_types"][1], 1)  # H

    def test_full_style_parsing(self) -> None:
        """Test parsing of full style LAMMPS data file with automatic detection."""
        self._test_style_parsing("full")

    def test_charge_style_parsing(self) -> None:
        """Test parsing of charge style LAMMPS data file with automatic detection."""
        self._test_style_parsing("charge")

    def test_bond_style_parsing(self) -> None:
        """Test parsing of bond style LAMMPS data file."""
        self._test_style_parsing("bond", explicit_style="bond")

    def test_full_style_no_comment_detection(self) -> None:
        """Test automatic detection of full style without style comment."""
        self._test_style_parsing("full_no_comment")

    def test_charge_style_no_comment_detection(self) -> None:
        """Test automatic detection of charge style without style comment."""
        self._test_style_parsing("charge_no_comment")

    def test_bond_style_no_comment_detection(self) -> None:
        """Test automatic detection of bond style without style comment."""
        self._test_style_parsing("bond_no_comment")

    def test_unsupported_atom_style(self) -> None:
        """Test that unsupported atom styles raise appropriate errors."""
        with self.assertRaises(ValueError) as context:
            dpdata.System(
                self.test_files["bond"],
                fmt="lammps/lmp",
                type_map=["O", "H"],
                atom_style="unsupported_style",
            )

        self.assertIn("Unsupported atom style", str(context.exception))

    def test_default_atomic_style(self) -> None:
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
