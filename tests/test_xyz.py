import unittest
import numpy as np
import tempfile
from context import dpdata
from comp_sys import CompSys, IsNoPBC

class TestToXYZ(unittest.TestCase):
    def test_to_xyz(self):
        with tempfile.NamedTemporaryFile('r') as f_xyz:
            dpdata.System(data={
                "atom_names": ["C", "O"],
                "atom_numbs": [1, 1],
                "atom_types": np.array([0, 1]),
                "coords": np.arange(6).reshape((1,2,3)),
                "cells": np.zeros((1,3,3)),
                "orig": np.zeros(3),
            }).to("xyz", f_xyz.name)
            xyz0 = f_xyz.read().strip()
        xyz1 = "2\n\nC 0.000000 1.000000 2.000000\nO 3.000000 4.000000 5.000000"
        self.assertEqual(xyz0, xyz1)


class TestFromXYZ(unittest.TestCase, CompSys, IsNoPBC):
    def setUp(self):
        self.places = 6
        # considering to_xyz has been tested..
        self.system_1 = dpdata.System(data={
                "atom_names": ["C", "O"],
                "atom_numbs": [1, 1],
                "atom_types": np.array([0, 1]),
                "coords": np.arange(6).reshape((1,2,3)),
                "cells": np.zeros((1,3,3)),
                "orig": np.zeros(3),
                "nopbc": True,
            })
        with tempfile.NamedTemporaryFile('r') as f_xyz:
            self.system_1.to("xyz", f_xyz.name)
            self.system_2 = dpdata.System(f_xyz.name, fmt="xyz")
