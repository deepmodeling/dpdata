from __future__ import annotations

import os
import pickle
import tempfile
import unittest

import numpy as np
from context import dpdata


class TestFeNNolFormat(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures with a simple water molecule system."""
        # Create a simple test system: water molecule (H2O)
        self.test_data = {
            "atom_names": ["H", "O"],
            "atom_numbs": [2, 1],
            "atom_types": np.array([0, 1, 0]),  # H, O, H
            "coords": np.array(
                [
                    [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]],  # frame 1
                    [[0.1, 0.0, 0.0], [0.0, 0.1, 1.0], [0.0, 1.1, 0.0]],  # frame 2
                ]
            ),  # 2 frames, 3 atoms, 3 coords
            "cells": np.array(
                [
                    [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],  # frame 1
                    [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],  # frame 2
                ]
            ),  # 2 frames, 3x3 cell
            "energies": np.array([-1.0, -1.1]),  # 2 frame energies in eV
            "forces": np.array(
                [
                    [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [-0.1, -0.1, 0.0]],  # frame 1
                    [[0.2, 0.0, 0.0], [0.0, 0.2, 0.0], [-0.2, -0.2, 0.0]],  # frame 2
                ]
            ),  # 2 frames, 3 atoms, 3 force components in eV/Angstrom
            "orig": np.array([0.0, 0.0, 0.0]),
            "nopbc": False,
        }

        self.system = dpdata.LabeledSystem(data=self.test_data)

    def test_fennol_export(self):
        """Test basic FeNNol format export functionality."""
        with tempfile.NamedTemporaryFile(suffix=".pkl", delete=False) as tmp_file:
            tmp_filename = tmp_file.name

        try:
            # Export to FeNNol format
            self.system.to("fennol", tmp_filename)

            # Check that file was created
            self.assertTrue(os.path.exists(tmp_filename))

            # Load and verify the FeNNol data
            with open(tmp_filename, "rb") as f:
                fennol_data = pickle.load(f)

            # Check main structure
            self.assertIn("training", fennol_data)
            self.assertIn("validation", fennol_data)
            self.assertIn("description", fennol_data)

            # Check that we have training and validation data
            training = fennol_data["training"]
            validation = fennol_data["validation"]

            # With default train_size=0.8 and 2 frames, we should have 1 training, 1 validation
            self.assertEqual(len(training), 1)
            self.assertEqual(len(validation), 1)

            # Check structure of training data
            sample = training[0]
            expected_keys = {
                "species",
                "coordinates",
                "formation_energy",
                "shifted_energy",
                "forces",
            }
            self.assertEqual(set(sample.keys()), expected_keys)

            # Check species
            expected_species = ["H", "O", "H"]
            self.assertEqual(sample["species"], expected_species)

            # Check coordinates (should be unchanged from Angstroms)
            np.testing.assert_array_almost_equal(
                sample["coordinates"], self.test_data["coords"][0]
            )

            # Check energy conversion (eV to kcal/mol)
            # 1 eV ≈ 23.06 kcal/mol
            expected_energy = self.test_data["energies"][0] * 23.06054783061903
            self.assertAlmostEqual(
                sample["formation_energy"], expected_energy, places=5
            )
            self.assertAlmostEqual(sample["shifted_energy"], expected_energy, places=5)

            # Check forces conversion
            expected_forces = self.test_data["forces"][0] * 23.06054783061903
            np.testing.assert_array_almost_equal(
                sample["forces"], expected_forces, decimal=5
            )

        finally:
            # Clean up
            if os.path.exists(tmp_filename):
                os.unlink(tmp_filename)

    def test_fennol_export_custom_train_size(self):
        """Test FeNNol export with custom training size."""
        with tempfile.NamedTemporaryFile(suffix=".pkl", delete=False) as tmp_file:
            tmp_filename = tmp_file.name

        try:
            # Export with train_size=0.5 (1 training, 1 validation from 2 frames)
            self.system.to("fennol", tmp_filename, train_size=0.5)

            with open(tmp_filename, "rb") as f:
                fennol_data = pickle.load(f)

            training = fennol_data["training"]
            validation = fennol_data["validation"]

            # Should have 1 training, 1 validation with train_size=0.5
            self.assertEqual(len(training), 1)
            self.assertEqual(len(validation), 1)

        finally:
            if os.path.exists(tmp_filename):
                os.unlink(tmp_filename)

    def test_fennol_export_all_training(self):
        """Test FeNNol export with all data as training."""
        with tempfile.NamedTemporaryFile(suffix=".pkl", delete=False) as tmp_file:
            tmp_filename = tmp_file.name

        try:
            # Export with train_size=1.0 (all training, no validation)
            self.system.to("fennol", tmp_filename, train_size=1.0)

            with open(tmp_filename, "rb") as f:
                fennol_data = pickle.load(f)

            training = fennol_data["training"]
            validation = fennol_data["validation"]

            # Should have 2 training, 0 validation
            self.assertEqual(len(training), 2)
            self.assertEqual(len(validation), 0)

        finally:
            if os.path.exists(tmp_filename):
                os.unlink(tmp_filename)

    def test_fennol_single_frame(self):
        """Test FeNNol export with single frame."""
        # Create single frame system
        single_frame_data = {
            k: v[:1]
            if k in ["coords", "cells", "energies"]
            else (v[:1] if k == "forces" else v)
            for k, v in self.test_data.items()
        }
        single_system = dpdata.LabeledSystem(data=single_frame_data)

        with tempfile.NamedTemporaryFile(suffix=".pkl", delete=False) as tmp_file:
            tmp_filename = tmp_file.name

        try:
            single_system.to("fennol", tmp_filename)

            with open(tmp_filename, "rb") as f:
                fennol_data = pickle.load(f)

            training = fennol_data["training"]
            validation = fennol_data["validation"]

            # With 1 frame and train_size=0.8, should have 0 training, 1 validation
            # (since int(1 * 0.8) = 0)
            self.assertEqual(len(training), 0)
            self.assertEqual(len(validation), 1)

        finally:
            if os.path.exists(tmp_filename):
                os.unlink(tmp_filename)


if __name__ == "__main__":
    unittest.main()
