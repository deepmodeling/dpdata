from __future__ import annotations

import unittest

import numpy as np
from context import dpdata


class TestMapAtomTypes(unittest.TestCase):
    def test_map_atom_types_preserves_current_atom_order(self):
        data = {
            "atom_names": ["O", "H"],
            "atom_numbs": [1, 2],
            "atom_types": np.array([1, 0, 1]),
            "orig": np.zeros(3),
            "cells": np.eye(3).reshape(1, 3, 3),
            "coords": np.zeros((1, 3, 3)),
        }

        system = dpdata.System(data=data)

        np.testing.assert_array_equal(
            system.map_atom_types({"H": 0, "O": 1}), np.array([0, 1, 0])
        )


class TestSetAtomTypes(unittest.TestCase):
    def setUp(self):
        self.system_1 = dpdata.LabeledSystem("poscars/vasprun.h2o.md.10.xml")
        self.type_1 = self.system_1.get_atom_types()
        self.system_types = np.array([0, 0, 1, 1, 1, 1])
        self.type_2 = self.system_1.map_atom_types(["H", "C", "O"])
        self.type_3 = self.system_1.map_atom_types({"H": 2, "C": 1, "O": 3})

    def test_types_func_1(self):
        atom_types = np.array([2, 2, 0, 0, 0, 0])
        atom_types_system_2 = self.type_2
        atom_types_system_1 = self.type_1
        for d0 in range(3):
            self.assertEqual(atom_types[d0], atom_types_system_2[d0])
        for d0 in range(3):
            self.assertEqual(self.system_types[d0], atom_types_system_1[d0])

    def test_types_func_2(self):
        atom_types = np.array([3, 3, 2, 2, 2, 2])
        atom_types_system_3 = self.type_3
        atom_types_system_1 = self.type_1
        for d0 in range(3):
            self.assertEqual(atom_types[d0], atom_types_system_3[d0])
        for d0 in range(3):
            self.assertEqual(self.system_types[d0], atom_types_system_1[d0])


if __name__ == "__main__":
    unittest.main()
