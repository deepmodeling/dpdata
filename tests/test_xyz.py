from __future__ import annotations

import tempfile
import unittest

import numpy as np
from comp_sys import CompSys, IsNoPBC
from context import dpdata


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


class TestExtXYZAlias(unittest.TestCase):
    def test_extxyz_alias_equivalence(self):
        """Test that extxyz alias produces same results as quip/gap/xyz."""
        # Use existing test file that works with quip/gap/xyz format
        test_file = "xyz/xyz_unittest.xyz"
        
        # Load with quip/gap/xyz format
        ms1 = dpdata.MultiSystems.from_file(test_file, "quip/gap/xyz")
        # Load with extxyz alias
        ms2 = dpdata.MultiSystems.from_file(test_file, "extxyz")
        
        # Should have same number of systems
        self.assertEqual(len(ms1.systems), len(ms2.systems))
        
        # System keys should match
        self.assertEqual(list(ms1.systems.keys()), list(ms2.systems.keys()))
        
        # Compare first system
        if ms1.systems:
            key = list(ms1.systems.keys())[0]
            sys1 = ms1.systems[key]
            sys2 = ms2.systems[key]
            
            # Basic properties should match
            self.assertEqual(len(sys1), len(sys2))
            self.assertEqual(sys1["atom_names"], sys2["atom_names"])
            np.testing.assert_array_equal(sys1["atom_types"], sys2["atom_types"])
            np.testing.assert_array_almost_equal(sys1["coords"], sys2["coords"])
            
            # Check that both have energies (labeled system property)
            if "energies" in sys1.data:
                np.testing.assert_array_almost_equal(sys1["energies"], sys2["energies"])
