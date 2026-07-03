from __future__ import annotations

import tempfile
import unittest

import numpy as np
from comp_sys import CompLabeledSys, CompSys, IsNoPBC
from context import dpdata

try:
    from ase.io import read, write
except ModuleNotFoundError:
    skip_ase = True
else:
    skip_ase = False


class TestToXYZ(unittest.TestCase):
    def test_to_xyz(self):
        with tempfile.NamedTemporaryFile("r") as f_xyz:
            dpdata.System(
                data={
                    "atom_names": ["C", "O"],
                    "atom_numbs": [1, 1],
                    "atom_types": np.array([0, 1]),
                    "coords": np.arange(6).reshape((1, 2, 3)),
                    "cells": np.zeros((1, 3, 3)),
                    "orig": np.zeros(3),
                }
            ).to("xyz", f_xyz.name)
            xyz0 = f_xyz.read().strip()
        xyz1 = "2\n\nC 0.000000 1.000000 2.000000\nO 3.000000 4.000000 5.000000"
        self.assertEqual(xyz0, xyz1)


class TestFromXYZ(unittest.TestCase, CompSys, IsNoPBC):
    def setUp(self):
        self.places = 6
        # considering to_xyz has been tested..
        self.system_1 = dpdata.System(
            data={
                "atom_names": ["C", "O"],
                "atom_numbs": [1, 1],
                "atom_types": np.array([0, 1]),
                "coords": np.arange(6).reshape((1, 2, 3)),
                "cells": np.zeros((1, 3, 3)),
                "orig": np.zeros(3),
                "nopbc": True,
            }
        )
        with tempfile.NamedTemporaryFile("r") as f_xyz:
            self.system_1.to("xyz", f_xyz.name)
            self.system_2 = dpdata.System(f_xyz.name, fmt="xyz")


@unittest.skipIf(skip_ase, "skip ASE related test. install ASE to fix")
class TestExtXYZASECrossCompatibility(unittest.TestCase):
    """Test cross-compatibility between dpdata extxyz and ASE extxyz."""

    def test_extxyz_format_compatibility_with_ase_read(self):
        """Test that dpdata's extxyz format can be read by ASE."""
        # Use existing test data that's known to work with dpdata extxyz parser
        test_file = "xyz/xyz_unittest.xyz"

        # First verify dpdata can read it
        multi_systems = dpdata.MultiSystems.from_file(test_file, fmt="extxyz")
        self.assertIsInstance(multi_systems, dpdata.MultiSystems)
        self.assertTrue(len(multi_systems.systems) > 0)

        # Test that ASE can also read the same file
        atoms_list = read(test_file, index=":", format="extxyz")
        self.assertIsInstance(atoms_list, list)
        self.assertTrue(len(atoms_list) > 0)

        # Check basic structure of first frame
        atoms = atoms_list[0]
        self.assertTrue(len(atoms) > 0)
        self.assertTrue(hasattr(atoms, "get_chemical_symbols"))

    def test_manual_extxyz_ase_to_dpdata(self):
        """Test cross-compatibility with a manually created compatible extxyz."""
        # Create a manually written extxyz content that should work with both
        extxyz_content = """2
energy=-10.5 Lattice="5.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 5.0" Properties=species:S:1:pos:R:3:Z:I:1:force:R:3
C 0.0 0.0 0.0 6 0.1 0.1 0.1
O 1.0 1.0 1.0 8 -0.1 -0.1 -0.1
"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".xyz", delete=False) as f:
            f.write(extxyz_content)
            f.flush()

            # Test with dpdata
            multi_systems = dpdata.MultiSystems.from_file(f.name, fmt="extxyz")
            self.assertIsInstance(multi_systems, dpdata.MultiSystems)
            self.assertTrue(len(multi_systems.systems) > 0)

            system_key = list(multi_systems.systems.keys())[0]
            system = multi_systems.systems[system_key]
            self.assertEqual(system.get_nframes(), 1)

            # Test with ASE (basic read)
            atoms = read(f.name, format="extxyz")
            self.assertEqual(len(atoms), 2)
            self.assertEqual(atoms.get_chemical_symbols(), ["C", "O"])

    def test_dpdata_xyz_to_ase_basic(self):
        """Test basic xyz reading between dpdata and ASE (simple compatibility check)."""
        # Create a simple xyz file using dpdata's basic xyz format
        simple_system = dpdata.System(
            data={
                "atom_names": ["C", "O"],
                "atom_numbs": [1, 1],
                "atom_types": np.array([0, 1]),
                "coords": np.array([[[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]]]),
                "cells": np.zeros((1, 3, 3)),
                "orig": np.zeros(3),
                "nopbc": True,
            }
        )

        with tempfile.NamedTemporaryFile(suffix=".xyz", mode="w+") as f:
            # Write basic xyz using dpdata
            simple_system.to("xyz", f.name)

            # Read with ASE
            atoms = read(f.name, format="xyz")

            # Verify basic structure
            self.assertEqual(len(atoms), 2)
            self.assertEqual(atoms.get_chemical_symbols(), ["C", "O"])

            # Check positions
            np.testing.assert_allclose(
                atoms.get_positions(), [[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]], rtol=1e-6
            )


@unittest.skipIf(skip_ase, "skip ASE related test. install ASE to fix")
class TestExtXYZEnergyForceCompatibility(unittest.TestCase, CompLabeledSys):
    """Test energy and force preservation between dpdata and ASE using CompLabeledSys."""

    def setUp(self):
        # Set precision for CompLabeledSys
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4

        # Create a manually written extxyz content with known energies and forces
        extxyz_content = """2
energy=-10.5 Lattice="5.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 5.0" Properties=species:S:1:pos:R:3:Z:I:1:force:R:3
C 0.0 1.0 2.0 6 0.1 0.1 0.1
O 3.0 4.0 5.0 8 -0.1 -0.1 -0.1
"""

        # Write the extxyz content to a file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".xyz", delete=False) as f:
            f.write(extxyz_content)
            f.flush()
            self.temp_file = f.name

        # Read with dpdata - this is our reference system
        multi_systems = dpdata.MultiSystems.from_file(self.temp_file, fmt="extxyz")
        system_key = list(multi_systems.systems.keys())[0]
        self.system_1 = multi_systems.systems[system_key]

        # Read with ASE
        atoms = read(self.temp_file, format="extxyz")

        # Write back to extxyz with ASE
        with tempfile.NamedTemporaryFile(suffix=".xyz", mode="w+", delete=False) as f2:
            self.temp_file2 = f2.name
            write(f2.name, atoms, format="extxyz")

        # Read back the ASE-written file with dpdata
        roundtrip_ms = dpdata.MultiSystems.from_file(self.temp_file2, fmt="extxyz")
        system_key = list(roundtrip_ms.systems.keys())[0]
        self.system_2 = roundtrip_ms.systems[system_key]

    def tearDown(self):
        import os

        try:
            os.unlink(self.temp_file)
            os.unlink(self.temp_file2)
        except (OSError, AttributeError):
            pass
