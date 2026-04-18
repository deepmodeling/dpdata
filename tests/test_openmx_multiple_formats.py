from __future__ import annotations

import unittest

import numpy as np
from context import dpdata


class TestOPENMXTRAJProps:
    def test_atom_names(self):
        self.assertEqual(self.system.data["atom_names"], ["Au", "H"])

    def test_atom_numbs(self):
        self.assertEqual(self.system.data["atom_numbs"], [3, 1])

    def test_atom_types(self):
        self.assertEqual(list(self.system.data["atom_types"]), [0, 0, 0, 1])

    def test_cell(self):
        cell = np.array(
            [[211.86272, 0.0, 0.0], [0.0, 2.88309, 0.0], [0.0, -1.44154, 2.49683]]
        )
        cells = np.array([cell, cell])
        self.assertEqual(self.system.get_nframes(), 2)
        for ff in range(self.system.get_nframes()):
            for ii in range(3):
                for jj in range(3):
                    self.assertAlmostEqual(
                        self.system["cells"][ff][ii][jj], cells[ff][ii][jj]
                    )

    def test_coord(self):
        coord = np.array(
            [
                [21.18627, 0.0, 1.66455],
                [23.5403, 1.44154, 0.83228],
                [25.89433, 0.0, 0.0],
                [28.24836, 0.0, 1.66455],
            ]
        )
        coords = np.array([coord, coord])
        for ff in range(self.system.get_nframes()):
            for ii in range(4):
                for jj in range(3):
                    self.assertAlmostEqual(
                        self.system["coords"][ff][ii][jj], coords[ff][ii][jj]
                    )


class TestOPENMXTraj(unittest.TestCase, TestOPENMXTRAJProps):
    def setUp(self):
        self.system = dpdata.System("openmx/Au111Surface", fmt="openmx/md")


class TestOPENMXLabeledTraj(unittest.TestCase, TestOPENMXTRAJProps):
    def setUp(self):
        self.system = dpdata.LabeledSystem("openmx/Au111Surface", fmt="openmx/md")


if __name__ == "__main__":
    unittest.main()
