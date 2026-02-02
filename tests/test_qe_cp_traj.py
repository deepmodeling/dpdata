from __future__ import annotations

import unittest

import numpy as np
from context import dpdata

bohr2ang = dpdata.unit.LengthConversion("bohr", "angstrom").value()


class TestCPTRAJProps:
    def test_atom_names(self):
        self.assertEqual(self.system.data["atom_names"], ["O", "H"])

    def test_atom_numbs(self):
        self.assertEqual(self.system.data["atom_numbs"], [64, 127])

    def test_atom_types(self):
        for ii in range(0, 64):
            self.assertEqual(self.system.data["atom_types"][ii], 0)
        for ii in range(64, 191):
            self.assertEqual(self.system.data["atom_types"][ii], 1)

    def test_cell(self):
        ref = bohr2ang * 23.5170 * np.eye(3)
        self.assertEqual(self.system.get_nframes(), 2)
        for ff in range(self.system.get_nframes()):
            for ii in range(3):
                for jj in range(3):
                    self.assertEqual(self.system["cells"][ff][ii][jj], ref[ii][jj])

    def test_coord(self):
        with open("qe.traj/oh-md.pos") as fp:
            lines = fp.read().rstrip("\n").split("\n")
        lines = lines[-191:]
        coords = []
        for ii in lines:
            coords.append([float(jj) for jj in ii.split()])
        coords = bohr2ang * np.array(coords)
        celll = bohr2ang * 23.5170
        for ii in range(coords.shape[0]):
            for jj in range(coords[ii].size):
                if coords[ii][jj] < 0:
                    coords[ii][jj] += celll
                elif coords[ii][jj] >= celll:
                    coords[ii][jj] -= celll
                self.assertAlmostEqual(
                    self.system["coords"][-1][ii][jj], coords[ii][jj]
                )


class TestCPTRAJTraj(unittest.TestCase, TestCPTRAJProps):
    def setUp(self):
        self.system = dpdata.System("qe.traj/oh-md", fmt="qe/cp/traj")


class TestCPTRAJLabeledTraj(unittest.TestCase, TestCPTRAJProps):
    def setUp(self):
        self.system = dpdata.LabeledSystem("qe.traj/oh-md", fmt="qe/cp/traj")


class TestConverCellDim(unittest.TestCase):
    def test_case_null(self):
        cell = dpdata.qe.traj.convert_celldm(8, [1, 1, 1])
        ref = np.eye(3)
        for ii in range(3):
            for jj in range(3):
                self.assertAlmostEqual(cell[ii][jj], ref[ii][jj])


class TestVirial(unittest.TestCase):
    def test(self):
        self.system = dpdata.LabeledSystem("qe.traj/si/si", fmt="qe/cp/traj")
        self.assertEqual(self.system["virials"].shape, (2, 3, 3))
        np.testing.assert_almost_equal(
            self.system["virials"][0],
            np.array(
                [
                    [0.31120718, -0.03261485, -0.02537362],
                    [-0.03261485, 0.3100397, 0.04211053],
                    [-0.02537362, 0.04211057, 0.30571264],
                ]
            ),
        )
        np.testing.assert_almost_equal(
            self.system["virials"][1],
            np.array(
                [
                    [0.31072979, -0.03151186, -0.02302297],
                    [-0.03151186, 0.30951293, 0.04078447],
                    [-0.02302297, 0.04078451, 0.30544987],
                ]
            ),
        )

    def test_raise(self):
        with self.assertRaises(RuntimeError) as c:
            self.system = dpdata.LabeledSystem(
                "qe.traj/si.wrongstr/si", fmt="qe/cp/traj"
            )
        self.assertTrue(
            "the step key between files are not consistent." in str(c.exception)
        )


if __name__ == "__main__":
    unittest.main()
