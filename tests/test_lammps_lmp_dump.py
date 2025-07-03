from __future__ import annotations

import os
import unittest

import numpy as np
from context import dpdata
from poscars.poscar_ref_oh import TestPOSCARoh

from dpdata.lammps.lmp import rotate_to_lower_triangle


class TestLmpDump(unittest.TestCase, TestPOSCARoh):
    def setUp(self):
        tmp_system = dpdata.System(
            os.path.join("poscars", "conf.lmp"), type_map=["O", "H"]
        )
        tmp_system.to_lammps_lmp("tmp.lmp")
        self.system = dpdata.System()
        self.system.from_lammps_lmp("tmp.lmp", type_map=["O", "H"])


class TestToFunc(unittest.TestCase, TestPOSCARoh):
    def setUp(self):
        tmp_system = dpdata.System(
            os.path.join("poscars", "conf.lmp"), type_map=["O", "H"]
        )
        tmp_system.to("lammps/lmp", "tmp.lmp")
        self.system = dpdata.System()
        self.system.from_fmt("tmp.lmp", fmt="lammps/lmp", type_map=["O", "H"])


class TestLmpRotateTriAngle(unittest.TestCase):
    def test_simple_cubic(self):
        cubic_cell = np.diag([5, 5, 5])
        cubic_coords = np.array([[0, 0, 0], [2.5, 2.5, 2.5]])
        rotated_cell, rotated_coords = rotate_to_lower_triangle(
            cubic_cell, cubic_coords
        )

        self.assertTrue(np.allclose(rotated_cell, np.diag([5, 5, 5])))
        self.assertTrue(np.allclose(rotated_coords, cubic_coords))

    def test_triclinic_cell(self):
        triclinic_cell = np.array([[4.0, 0.0, 0.0], [1.0, 3.0, 0.0], [0.5, 0.5, 2.5]])

        triclinic_coords = np.array([[0.5, 0.5, 0.5], [3.0, 2.0, 1.0]])
        rotated_cell, rotated_coords = rotate_to_lower_triangle(
            triclinic_cell.copy(), triclinic_coords.copy()
        )

        self.assertTrue(np.isclose(rotated_cell[0, 1], 0, atol=1e-10))
        self.assertTrue(np.isclose(rotated_cell[0, 2], 0, atol=1e-10))
        self.assertTrue(np.isclose(rotated_cell[1, 2], 0, atol=1e-10))

        self.assertTrue(rotated_cell[0, 0] > 0)
        self.assertTrue(rotated_cell[1, 1] > 0)
        self.assertTrue(rotated_cell[2, 2] > 0)

        # check the distance between two atoms
        self.assertTrue(
            np.isclose(
                np.linalg.norm(rotated_coords[0] - rotated_coords[1]),
                np.linalg.norm(triclinic_coords[0] - triclinic_coords[1]),
                atol=1e-10,
            )
        )

    def test_negative_diagonal(self):
        negative_cell = np.array([[-3.0, 1.0, 0.5], [0.0, -2.0, 0.3], [0.0, 0.0, -4.0]])
        negative_coords = np.array([[1.0, 0.5, -1.0], [0.5, 1.0, -2.0]])
        rotated_cell, rotated_coords = rotate_to_lower_triangle(
            negative_cell.copy(), negative_coords.copy()
        )

        self.assertTrue(np.isclose(rotated_cell[0, 1], 0, atol=1e-10))
        self.assertTrue(np.isclose(rotated_cell[0, 2], 0, atol=1e-10))
        self.assertTrue(np.isclose(rotated_cell[1, 2], 0, atol=1e-10))

        self.assertTrue(rotated_cell[0, 0] > 0)
        self.assertTrue(rotated_cell[1, 1] > 0)
        self.assertTrue(rotated_cell[2, 2] > 0)

        self.assertTrue(
            np.isclose(
                np.linalg.norm(negative_coords[0] - negative_coords[1]),
                np.linalg.norm(rotated_coords[0] - rotated_coords[1]),
                atol=1e-6,
            )
        )


if __name__ == "__main__":
    unittest.main()
