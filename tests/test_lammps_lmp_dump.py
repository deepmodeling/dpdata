from __future__ import annotations

import os
import tempfile
import unittest

import numpy as np
from context import dpdata
from poscars.poscar_ref_oh import TestPOSCARoh

from dpdata.lammps.lmp import rotate_to_lower_triangle

TEST_DIR = os.path.dirname(__file__)
POSCAR_CONF_LMP = os.path.join(TEST_DIR, "poscars", "conf.lmp")


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


class TestDumpToFunc(unittest.TestCase, TestPOSCARoh):
    def setUp(self):
        tmp_system = dpdata.System(
            os.path.join("poscars", "conf.lmp"), type_map=["O", "H"]
        )
        tmp_system.to("lammps/dump", "tmp.dump")
        self.system = dpdata.System()
        self.system.from_fmt("tmp.dump", fmt="lammps/dump", type_map=["O", "H"])


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


class TestLmpDumpMasses(unittest.TestCase):
    def test_dump_known_elements_writes_masses(self):
        system = dpdata.System(POSCAR_CONF_LMP, type_map=["O", "H"])
        with tempfile.TemporaryDirectory() as tmpdir:
            output = os.path.join(tmpdir, "tmp_masses.lmp")
            system.to_lammps_lmp(output)
            with open(output) as f:
                content = f.read()

        self.assertIn("Masses\n", content)
        self.assertIn("     1   15.9994000000 # O", content)
        self.assertIn("     2    1.0079400000 # H", content)
        self.assertLess(content.index("Masses\n"), content.index("Atoms # atomic\n"))

    def test_dump_unknown_types_skips_masses(self):
        system = dpdata.System(POSCAR_CONF_LMP)
        with tempfile.TemporaryDirectory() as tmpdir:
            output = os.path.join(tmpdir, "tmp_unknown_types.lmp")
            system.to_lammps_lmp(output)
            with open(output) as f:
                content = f.read()

        self.assertNotIn("Masses\n", content)

    def test_dump_rejects_mismatched_explicit_masses(self):
        system = dpdata.System(POSCAR_CONF_LMP, type_map=["O", "H"])
        system.data["masses"] = np.array([15.9994, 1.00794, 99.0])

        with tempfile.TemporaryDirectory() as tmpdir:
            output = os.path.join(tmpdir, "tmp_bad_masses.lmp")
            with self.assertRaisesRegex(ValueError, r'system\["masses"\]'):
                system.to_lammps_lmp(output)


if __name__ == "__main__":
    unittest.main()
