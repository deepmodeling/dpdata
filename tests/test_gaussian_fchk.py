from __future__ import annotations

import unittest

import numpy as np
from context import dpdata


class TestGaussianLoadFchk(unittest.TestCase):
    def setUp(self):
        self.system = dpdata.LabeledSystem(
            "gaussian/waterfreq.gaussianfchk", fmt="gaussian/fchk"
        )
        self.atom_names = ["H", "O"]
        self.atom_numbs = [2, 1]
        self.nframes = 1
        self.atom_types = [1, 0, 0]

    def test_atom_names(self):
        """Test that atom names are correctly read."""
        self.assertEqual(self.system.data["atom_names"], self.atom_names)

    def test_atom_numbs(self):
        """Test that atom numbers are correctly read."""
        self.assertEqual(self.system.data["atom_numbs"], self.atom_numbs)

    def test_nframes(self):
        """Test that number of frames is correct."""
        self.assertEqual(len(self.system), self.nframes)

    def test_atom_types(self):
        """Test that atom types are correctly assigned."""
        for ii in range(len(self.atom_types)):
            self.assertEqual(self.system.data["atom_types"][ii], self.atom_types[ii])

    def test_nopbc(self):
        """Test that nopbc is set to True for fchk files."""
        self.assertEqual(self.system.nopbc, True)

    def test_energies_exist(self):
        """Test that energies are present and have correct shape."""
        self.assertIn("energies", self.system.data)
        self.assertEqual(len(self.system.data["energies"]), 1)
        self.assertIsInstance(self.system.data["energies"], np.ndarray)

    def test_forces_exist(self):
        """Test that forces are present and have correct shape."""
        self.assertIn("forces", self.system.data)
        self.assertEqual(self.system.data["forces"].shape, (1, 3, 3))
        self.assertIsInstance(self.system.data["forces"], np.ndarray)

    def test_hessian_exist(self):
        """Test that hessian matrix is present and has correct shape."""
        self.assertIn("hessian", self.system.data)
        self.assertEqual(self.system.data["hessian"].shape, (1, 3, 3, 3, 3))
        self.assertIsInstance(self.system.data["hessian"], np.ndarray)

    def test_coords_exist(self):
        """Test that coordinates are present and have correct shape."""
        self.assertIn("coords", self.system.data)
        self.assertEqual(self.system.data["coords"].shape, (1, 3, 3))
        self.assertIsInstance(self.system.data["coords"], np.ndarray)

    def test_cells_exist(self):
        """Test that cells are present and have correct shape."""
        self.assertIn("cells", self.system.data)
        self.assertEqual(self.system.data["cells"].shape, (1, 3, 3))
        self.assertIsInstance(self.system.data["cells"], np.ndarray)

    def test_orig_exist(self):
        """Test that origin is present and has correct shape."""
        self.assertIn("orig", self.system.data)
        self.assertEqual(self.system.data["orig"].shape, (3,))
        self.assertIsInstance(self.system.data["orig"], np.ndarray)


class TestGaussianFchkVsLog(unittest.TestCase):
    """Test to compare results from fchk and log files."""

    def setUp(self):
        # Load both fchk and log files
        self.system_fchk = dpdata.LabeledSystem(
            "gaussian/waterfreq.gaussianfchk", fmt="gaussian/fchk"
        )
        self.system_log = dpdata.LabeledSystem(
            "gaussian/waterfreq.gaussianlog", fmt="gaussian/log"
        )

        # Get conversion factors from dpdata
        from dpdata.unit import EnergyConversion, ForceConversion, LengthConversion

        self.energy_convert = EnergyConversion("hartree", "eV").value()
        self.force_convert = ForceConversion("hartree/bohr", "eV/angstrom").value()
        self.length_convert = LengthConversion("bohr", "angstrom").value()

    def test_energies_consistency(self):
        """Test that energies from fchk and log files are consistent."""
        # Check that both files have energies
        self.assertIn("energies", self.system_fchk.data)
        self.assertIn("energies", self.system_log.data)

        # Check that energies have the same length
        self.assertEqual(
            len(self.system_fchk.data["energies"]),
            len(self.system_log.data["energies"]),
        )

        # Check that energies are equal (allowing for small numerical differences)
        fchk_energy = self.system_fchk.data["energies"][0]
        log_energy = self.system_log.data["energies"][0]
        self.assertAlmostEqual(fchk_energy, log_energy, places=6)

    def test_forces_consistency(self):
        """Test that forces from fchk and log files are consistent."""
        # Check that both files have forces
        self.assertIn("forces", self.system_fchk.data)
        self.assertIn("forces", self.system_log.data)

        # Check that forces have the same shape
        self.assertEqual(
            self.system_fchk.data["forces"].shape, self.system_log.data["forces"].shape
        )

        # Check that forces are equal (allowing for small numerical differences)
        fchk_forces = self.system_fchk.data["forces"][0]
        log_forces = self.system_log.data["forces"][0]
        np.testing.assert_array_almost_equal(fchk_forces, log_forces, decimal=6)

    def test_coordinates_consistency(self):
        """Test that coordinates from fchk and log files are consistent."""
        # Check that both files have coordinates
        self.assertIn("coords", self.system_fchk.data)
        self.assertIn("coords", self.system_log.data)

        # Check that coordinates have the same shape
        self.assertEqual(
            self.system_fchk.data["coords"].shape, self.system_log.data["coords"].shape
        )

        # Check that coordinates are equal (allowing for small numerical differences)
        fchk_coords = self.system_fchk.data["coords"][0]
        log_coords = self.system_log.data["coords"][0]
        np.testing.assert_array_almost_equal(fchk_coords, log_coords, decimal=6)

    def test_atom_info_consistency(self):
        """Test that atom information is consistent between fchk and log files."""
        # Check atom names
        self.assertEqual(
            self.system_fchk.data["atom_names"], self.system_log.data["atom_names"]
        )

        # Check atom numbers
        self.assertEqual(
            self.system_fchk.data["atom_numbs"], self.system_log.data["atom_numbs"]
        )

        # Check atom types
        np.testing.assert_array_equal(
            self.system_fchk.data["atom_types"], self.system_log.data["atom_types"]
        )

    def test_system_properties_consistency(self):
        """Test that system properties are consistent between fchk and log files."""
        # Check number of frames
        self.assertEqual(len(self.system_fchk), len(self.system_log))

        # Check nopbc property
        self.assertEqual(self.system_fchk.nopbc, self.system_log.nopbc)

        # Check that both have the same data keys
        fchk_keys = set(self.system_fchk.data.keys())
        log_keys = set(self.system_log.data.keys())

        # fchk has hessian, log doesn't, so we exclude it from comparison
        fchk_keys.discard("hessian")

        self.assertEqual(fchk_keys, log_keys)


if __name__ == "__main__":
    unittest.main()
