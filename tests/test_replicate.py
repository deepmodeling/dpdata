from __future__ import annotations

import unittest

import numpy as np
from comp_sys import CompSys, IsPBC
from context import dpdata


class TestReplicate123(unittest.TestCase, CompSys, IsPBC):
    def setUp(self):
        system_1_origin = dpdata.System("poscars/POSCAR.SiC", fmt="vasp/poscar")
        self.system_1 = system_1_origin.replicate(
            (
                1,
                2,
                3,
            )
        )
        self.system_2 = dpdata.System(
            "poscars/POSCAR.SiC.replicate123", fmt="vasp/poscar"
        )
        self.places = 6


class TestReplicate123_not_change_origin(unittest.TestCase, CompSys, IsPBC):
    def setUp(self):
        self.system_1 = dpdata.System("poscars/POSCAR.SiC", fmt="vasp/poscar")
        self.system_1.replicate(
            (
                1,
                2,
                3,
            )
        )
        self.system_2 = dpdata.System("poscars/POSCAR.SiC", fmt="vasp/poscar")
        self.places = 6

class TestReplicateTriclinicBox(unittest.TestCase, CompSys, IsPBC):
    def setUp(self):
        self.system_1 = dpdata.System()
        self.system_1.data["atom_names"] = ["foo", "bar"]
        self.system_1.data["atom_types"] = np.array([1, 0], dtype=int)
        self.system_1.data["atom_numbs"] = [1, 1]
        self.system_1.data["cells"] = np.array([10, 0, 0, 0, 10, 0, 0, 0, 10], dtype=float).reshape(1,3,3)
        self.system_1.data["coords"] = np.array([0, 0, 0, 0, 0, 1], dtype=float).reshape(1,2,3)
        self.system_1 = self.system_1.replicate([2,1,1])

        self.system_2 = dpdata.System()
        self.system_2.data["atom_names"] = ["foo", "bar"]
        self.system_2.data["atom_types"] = np.array([1, 1, 0, 0], dtype=int)
        self.system_2.data["atom_numbs"] = [2, 2]
        self.system_2.data["cells"] = np.array([20, 0, 0, 0, 10, 0, 0, 0, 10], dtype=float).reshape(1,3,3)
        self.system_2.data["coords"] = np.array([0, 0, 0, 10, 0, 0, 0, 0, 1, 10, 0, 1], dtype=float).reshape(1,4,3)
        self.places = 6


if __name__ == "__main__":
    unittest.main()
