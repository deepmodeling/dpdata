from __future__ import annotations

import unittest

import numpy as np
from context import dpdata

from dpdata.formats.qe.traj import convert_celldm, load_cell_parameters

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


class TestCPTRAJInputCellUnits(unittest.TestCase):
    def test_angstrom_cell_without_cel_trajectory(self):
        # Earlier tests did not expose this regression: the only missing-.cel
        # fixture used ibrav/celldm, whose cell is correctly expressed in Bohr.
        # Exercise the distinct fallback path where CELL_PARAMETERS already
        # stores angstrom values and therefore must not be converted again.
        system = dpdata.System(
            "qe.traj/angstrom_no_cel/cp",
            fmt="qe/cp/traj",
        )

        np.testing.assert_allclose(
            system["cells"][0],
            np.eye(3) * 19.7299995422,
        )

    def test_bohr_cell_without_cel_trajectory(self):
        system = dpdata.System(
            "qe.traj/bohr_no_cel/cp",
            fmt="qe/cp/traj",
        )

        np.testing.assert_allclose(
            system["cells"][0],
            np.eye(3) * 2.0 * bohr2ang,
        )

    def test_alat_cell_without_cel_trajectory(self):
        system = dpdata.System(
            "qe.traj/alat_no_cel/cp",
            fmt="qe/cp/traj",
        )

        np.testing.assert_allclose(
            system["cells"][0],
            np.eye(3) * 10.0 * bohr2ang,
        )

    def test_omitted_unit_uses_a_lattice_parameter(self):
        system = dpdata.System(
            "qe.traj/a_no_cel/cp",
            fmt="qe/cp/traj",
        )
        np.testing.assert_allclose(system["cells"][0], np.eye(3) * 5.5)

    def test_alat_without_lattice_parameter_raises(self):
        with self.assertRaisesRegex(ValueError, "requires celldm\\(1\\) or A"):
            load_cell_parameters(["CELL_PARAMETERS {alat}", "1 0 0", "0 1 0", "0 0 1"])

    def test_omitted_unit_without_lattice_parameter_raises(self):
        with self.assertRaisesRegex(ValueError, "without a unit requires"):
            load_cell_parameters(["CELL_PARAMETERS", "1 0 0", "0 1 0", "0 0 1"])

    def test_unsupported_cell_unit_raises(self):
        with self.assertRaisesRegex(ValueError, "unsupported CELL_PARAMETERS unit"):
            load_cell_parameters(
                ["CELL_PARAMETERS {crystal}", "1 0 0", "0 1 0", "0 0 1"]
            )

    def test_ambiguous_cell_unit_raises(self):
        with self.assertRaisesRegex(ValueError, "ambiguous CELL_PARAMETERS unit"):
            load_cell_parameters(
                ["CELL_PARAMETERS {alat} bohr", "1 0 0", "0 1 0", "0 0 1"],
                lattice_parameter=5.0,
            )

    def test_missing_cell_parameters_for_ibrav_zero_raises(self):
        with self.assertRaisesRegex(
            ValueError, "CELL_PARAMETERS is required when ibrav is 0"
        ):
            dpdata.System(
                "qe.traj/missing_cell_no_cel/cp",
                fmt="qe/cp/traj",
            )


class TestConverCellDim(unittest.TestCase):
    # Reference cells follow latgen in QE's Modules/latgen.f90. They were
    # checked against each lattice family's lengths, angles, and volume.
    CELLDIM = np.array([4.0, 1.5, 2.0, 0.5, 0.4, 0.6])
    GOLDEN = {
        1: [[4, 0, 0], [0, 4, 0], [0, 0, 4]],
        2: [[-2, 0, 2], [0, 2, 2], [-2, 2, 0]],
        3: [[2, 2, 2], [-2, 2, 2], [-2, -2, 2]],
        -3: [[-2, 2, 2], [2, -2, 2], [2, 2, -2]],
        4: [[4, 0, 0], [-2, 2 * np.sqrt(3), 0], [0, 0, 8]],
        5: [
            [2, -2 / np.sqrt(3), 4 * np.sqrt(2 / 3)],
            [0, 4 / np.sqrt(3), 4 * np.sqrt(2 / 3)],
            [-2, -2 / np.sqrt(3), 4 * np.sqrt(2 / 3)],
        ],
        -5: [
            [0, 2 * np.sqrt(2), 2 * np.sqrt(2)],
            [2 * np.sqrt(2), 0, 2 * np.sqrt(2)],
            [2 * np.sqrt(2), 2 * np.sqrt(2), 0],
        ],
        6: [[4, 0, 0], [0, 4, 0], [0, 0, 8]],
        7: [[2, -2, 4], [2, 2, 4], [-2, -2, 4]],
        8: [[4, 0, 0], [0, 6, 0], [0, 0, 8]],
        9: [[2, 3, 0], [-2, 3, 0], [0, 0, 8]],
        -9: [[2, -3, 0], [2, 3, 0], [0, 0, 8]],
        91: [[4, 0, 0], [0, 3, -4], [0, 3, 4]],
        10: [[2, 0, 4], [2, 3, 0], [0, 3, 4]],
        11: [[2, 3, 4], [-2, 3, 4], [-2, -3, 4]],
        12: [[4, 0, 0], [3, 3 * np.sqrt(3), 0], [0, 0, 8]],
        -12: [[4, 0, 0], [0, 6, 0], [3.2, 0, 8 * np.sqrt(1 - 0.4**2)]],
        13: [[2, 0, -4], [3, 3 * np.sqrt(3), 0], [2, 0, 4]],
        -13: [[2, 3, 0], [-2, 3, 0], [3.2, 0, 8 * np.sqrt(1 - 0.4**2)]],
        14: [
            [4, 0, 0],
            [3.6, 4.8, 0],
            [
                3.2,
                2.6,
                8
                * np.sqrt(
                    (1.0 + 2.0 * 0.5 * 0.4 * 0.6 - 0.5**2 - 0.4**2 - 0.6**2)
                    / (1.0 - 0.6**2)
                ),
            ],
        ],
    }

    def test_all_ibrav(self):
        for ibrav, expected in self.GOLDEN.items():
            with self.subTest(ibrav=ibrav):
                cell = convert_celldm(ibrav, self.CELLDIM)
                np.testing.assert_allclose(cell, expected, atol=1e-7)

    def test_lengths_angles_and_volume(self):
        # (lens, sorted cos of the three pairs, |det|), independent of axis choice.
        expected = {
            1: ([4, 4, 4], [0, 0, 0], 64),
            2: ([4 / np.sqrt(2)] * 3, [0.5, 0.5, 0.5], 16),
            4: ([4, 4, 8], [-0.5, 0, 0], 0.5 * np.sqrt(3) * 4 * 4 * 8),
            5: ([4, 4, 4], [0.5, 0.5, 0.5], 45.25483),
            8: ([4, 6, 8], [0, 0, 0], 192),
            12: ([4, 6, 8], [0, 0, 0.5], 192 * np.sqrt(0.75)),
            -12: ([4, 6, 8], [0, 0.4, 0], 192 * np.sqrt(1 - 0.4**2)),
            14: ([4, 6, 8], [0.4, 0.5, 0.6], 131.6286),
        }
        for ibrav, (lens, cos, vol) in expected.items():
            with self.subTest(ibrav=ibrav):
                cell = convert_celldm(ibrav, self.CELLDIM)
                np.testing.assert_allclose(
                    np.linalg.norm(cell, axis=1), lens, atol=1e-6
                )
                pairs = [
                    np.dot(cell[i], cell[j])
                    / np.linalg.norm(cell[i])
                    / np.linalg.norm(cell[j])
                    for i in range(3)
                    for j in range(i + 1, 3)
                ]
                np.testing.assert_allclose(sorted(pairs), sorted(cos), atol=1e-6)
                self.assertAlmostEqual(abs(np.linalg.det(cell)), vol, places=4)

    def test_unsupported_ibrav_raises(self):
        for ibrav in (0, 15, -1, 99):
            with self.subTest(ibrav=ibrav):
                with self.assertRaises(RuntimeError):
                    convert_celldm(ibrav, self.CELLDIM)

    def test_invalid_celldm_raises(self):
        # (ibrav, celldm index to invalidate, replacement, QE-slot label)
        bad_inputs = [
            (1, 0, 0.0, "celldm(1)=0"),
            (1, 0, -1.0, "celldm(1)<0"),
            (4, 2, 0.0, "ibrav=4 celldm(3)<=0"),
            (4, 2, -1.0, "ibrav=4 celldm(3)<0"),
            (5, 3, -0.5, "ibrav=5 cosAB=-0.5"),
            (5, 3, 1.0, "ibrav=5 cosAB=1"),
            (-5, 3, 1.5, "ibrav=-5 cosAB out of range"),
            (6, 2, 0.0, "ibrav=6 celldm(3)<=0"),
            (7, 2, -0.1, "ibrav=7 celldm(3)<0"),
            (8, 1, 0.0, "ibrav=8 celldm(2)=0"),
            (8, 2, -1.0, "ibrav=8 celldm(3)<0"),
            (9, 1, 0.0, "ibrav=9 celldm(2)=0"),
            (-9, 2, 0.0, "ibrav=-9 celldm(3)=0"),
            (91, 1, 0.0, "ibrav=91 celldm(2)=0"),
            (10, 2, 0.0, "ibrav=10 celldm(3)=0"),
            (11, 1, 0.0, "ibrav=11 celldm(2)=0"),
            (12, 3, 1.0, "ibrav=12 cosAB=+1 (degenerate sen)"),
            (12, 3, -1.0, "ibrav=12 cosAB=-1"),
            (12, 1, 0.0, "ibrav=12 celldm(2)=0"),
            (-12, 4, 1.0, "ibrav=-12 cosAC=+1"),
            (-12, 2, -1.0, "ibrav=-12 celldm(3)<0"),
            (13, 3, 1.0, "ibrav=13 cosAB=+1"),
            (-13, 4, -1.0, "ibrav=-13 cosAC=-1"),
            (14, 1, 0.0, "ibrav=14 celldm(2)=0"),
            (14, 2, 0.0, "ibrav=14 celldm(3)=0"),
            (14, 3, 1.0, "ibrav=14 cosBC=+1"),
            (14, 4, -1.0, "ibrav=14 cosAC=-1"),
            (14, 5, 1.0, "ibrav=14 cosAB=+1 (singam=0)"),
        ]
        for ibrav, idx, value, label in bad_inputs:
            with self.subTest(case=label):
                celldm = self.CELLDIM.copy()
                celldm[idx] = value
                with self.assertRaises(RuntimeError):
                    convert_celldm(ibrav, celldm)

    def test_ibrav_14_degenerate_term_rejects_zero(self):
        # cosBC=cosAC=cosAB=0.5 gives 1 + 2*(0.5)**3 - 3*(0.5)**2 = 0.875 > 0 (valid)
        # But cosAB=0, cosAC=1-eps, cosBC=1-eps pushes term → 0 exactly at boundary
        celldm = self.CELLDIM.copy()
        celldm[3] = 0.0  # cosBC
        celldm[4] = 1.0 - 1e-9  # cosAC (just inside range)
        celldm[5] = 1.0 - 1e-9  # cosAB (just inside range)
        # term = 1 + 0 - 0 - (1-1e-9)^2 - (1-1e-9)^2 ≈ -1 + 4e-9 → still rejected
        with self.assertRaises(RuntimeError):
            convert_celldm(14, celldm)


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
