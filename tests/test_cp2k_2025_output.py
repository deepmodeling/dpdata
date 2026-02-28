from __future__ import annotations

import unittest

import numpy as np
from comp_sys import CompLabeledSys
from context import dpdata


class TestCp2k2025Output(unittest.TestCase, CompLabeledSys):
    """Test CP2K 2025 output format parsing."""

    def setUp(self):
        self.system_1 = dpdata.LabeledSystem(
            "cp2k/cp2k_2025_output/cp2k_2025_output", fmt="cp2k/output"
        )
        self.system_2 = dpdata.LabeledSystem(
            "cp2k/cp2k_2025_output/deepmd", fmt="deepmd/npy"
        )
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4

    def test_energy_extraction(self):
        """Test that energy is correctly extracted from CP2K 2025 format."""
        # Energy should be -7.364190264587725 hartree
        # Using dpdata's conversion factor: -200.3898256786414 eV
        expected_energy = -200.3898256786414
        self.assertAlmostEqual(
            self.system_1.data["energies"][0],
            expected_energy,
            places=5
        )

    def test_forces_extraction(self):
        """Test that forces are correctly extracted from CP2K 2025 format."""
        # Forces should be converted from hartree/bohr to eV/angstrom
        self.assertEqual(self.system_1.data["forces"].shape, (1, 2, 3))
        # Check first atom force x-component
        self.assertAlmostEqual(
            self.system_1.data["forces"][0][0][0],
            -2.94874881,
            places=5
        )


class TestCp2k2023FormatStillWorks(unittest.TestCase, CompLabeledSys):
    """Test that CP2K 2023 format still works (regression test)."""

    def setUp(self):
        self.system_1 = dpdata.LabeledSystem(
            "cp2k/cp2k_normal_output/cp2k_output", fmt="cp2k/output"
        )
        self.system_2 = dpdata.LabeledSystem(
            "cp2k/cp2k_normal_output/deepmd", fmt="deepmd/npy"
        )
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4


if __name__ == "__main__":
    unittest.main()
