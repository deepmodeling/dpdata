from __future__ import annotations

import os
import unittest

import numpy as np
from context import dpdata

from dpdata.unit import EnergyConversion, ForceConversion, LengthConversion

length_convert = LengthConversion("bohr", "angstrom").value()
energy_convert = EnergyConversion("hartree", "eV").value()
force_convert = ForceConversion("hartree/bohr", "eV/angstrom").value()


class TestN2P2(unittest.TestCase):
    def setUp(self):
        self.data_ref = {
            "atom_numbs": [1, 2],
            "atom_names": ["O", "H"],
            "atom_types": np.array([0, 1, 1]),
            "orig": np.array([0.0, 0.0, 0.0]),
            "cells": np.array(
                [
                    [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
                    [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
                ]
            ),
            "coords": np.array(
                [
                    [[1.0, 0.0, 0.0], [0.0, 0.0, 1.5], [1.0, 0.0, 3.0]],
                    [[2.0, 1.0, 1.0], [1.0, 1.0, 2.5], [2.0, 1.0, 4.0]],
                ]
            ),
            "energies": np.array([1.2, 2.3]),
            "forces": np.array(
                [
                    [[0.5, 0.0, 0.0], [0.0, 0.0, 0.75], [0.5, 0.0, 1.5]],
                    [[2.5, 2.0, 2.0], [2.0, 2.0, 2.75], [2.5, 2.0, 3.5]],
                ]
            ),
        }

    def test_n2p2_from_labeled_system(self):
        data = dpdata.LabeledSystem("n2p2/input.data", fmt="n2p2")
        for key in self.data_ref:
            if key == "atom_numbs":
                self.assertEqual(data[key], self.data_ref[key])
            elif key == "atom_names":
                self.assertEqual(data[key], self.data_ref[key])
            elif key == "atom_types":
                np.testing.assert_array_equal(data[key], self.data_ref[key])
            else:
                np.testing.assert_array_almost_equal(
                    data[key], self.data_ref[key], decimal=5
                )

    def test_n2p2_to_labeled_system(self):
        output_file = "n2p2/output.data"
        data = dpdata.LabeledSystem.from_dict({"data": self.data_ref})
        data.to_n2p2(output_file)
        ref_file = "n2p2/input.data"
        with open(ref_file) as file1, open(output_file) as file2:
            file1_lines = file1.readlines()
            file2_lines = file2.readlines()

        file1_lines = [line.strip("\n") for line in file1_lines]
        file2_lines = [line.strip("\n") for line in file2_lines]

        self.assertListEqual(file1_lines, file2_lines)

    def tearDown(self):
        if os.path.isfile("n2p2/output.data"):
            os.remove("n2p2/output.data")
