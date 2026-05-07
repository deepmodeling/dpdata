from __future__ import annotations

import os
import shutil
import tempfile
import unittest

import numpy as np
from comp_sys import CompLabeledSys, IsPBC
from context import dpdata

from dpdata.formats.amber.md import cell_lengths_angles_to_cell

try:
    import parmed  # noqa: F401
except ModuleNotFoundError:
    skip_parmed_related_test = True
else:
    skip_parmed_related_test = False


class TestAmberMD(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        self.system_1 = dpdata.LabeledSystem("amber/02_Heat", fmt="amber/md")
        self.system_1.to("deepmd/npy", "tmp.deepmd.npy")
        self.system_2 = dpdata.LabeledSystem("tmp.deepmd.npy", fmt="deepmd/npy")
        self.places = 5
        self.e_places = 4
        self.f_places = 6
        self.v_places = 6

    def tearDown(self):
        if os.path.exists("tmp.deepmd.npy"):
            shutil.rmtree("tmp.deepmd.npy")


class TestAmberMDNonOrthogonal(unittest.TestCase):
    def test_cell_lengths_angles_to_cell(self):
        cell_lengths = np.array([[10.0, 10.0, 15.0], [8.0, 10.0, 12.0]])
        cell_angles = np.array([[90.0, 90.0, 120.0], [70.0, 80.0, 110.0]])

        cells = cell_lengths_angles_to_cell(cell_lengths, cell_angles)

        self.assertEqual(cells.shape, (2, 3, 3))
        np.testing.assert_allclose(
            cells[0],
            np.array(
                [
                    [10.0, 0.0, 0.0],
                    [-5.0, 5.0 * np.sqrt(3.0), 0.0],
                    [0.0, 0.0, 15.0],
                ]
            ),
            atol=1e-12,
        )
        for frame, lengths, angles in zip(cells, cell_lengths, cell_angles):
            np.testing.assert_allclose(np.linalg.norm(frame, axis=1), lengths)
            alpha = np.rad2deg(
                np.arccos(
                    np.dot(frame[1], frame[2])
                    / (np.linalg.norm(frame[1]) * np.linalg.norm(frame[2]))
                )
            )
            beta = np.rad2deg(
                np.arccos(
                    np.dot(frame[0], frame[2])
                    / (np.linalg.norm(frame[0]) * np.linalg.norm(frame[2]))
                )
            )
            gamma = np.rad2deg(
                np.arccos(
                    np.dot(frame[0], frame[1])
                    / (np.linalg.norm(frame[0]) * np.linalg.norm(frame[1]))
                )
            )
            np.testing.assert_allclose([alpha, beta, gamma], angles)

    def test_invalid_cell_lengths(self):
        cell_lengths = np.array([[0.0, 8.0, 12.0], [5.0, -8.0, 12.0]])
        cell_angles = np.array([[90.0, 90.0, 90.0], [90.0, 90.0, 90.0]])

        with self.assertRaisesRegex(RuntimeError, "Invalid AMBER cell lengths"):
            cell_lengths_angles_to_cell(cell_lengths, cell_angles)

    def test_invalid_cell_angles(self):
        cell_lengths = np.array(
            [
                [5.0, 8.0, 12.0],
                [5.0, 8.0, 12.0],
                [5.0, 8.0, 12.0],
                [5.0, 8.0, 12.0],
            ]
        )
        cell_angles = np.array(
            [
                [60.0, 70.0, 130.0],
                [90.0, 90.0, 0.0],
                [90.0, 90.0, 180.0],
                [90.0, 90.0, 180.1],
            ]
        )

        with self.assertRaisesRegex(RuntimeError, "Invalid AMBER cell angles"):
            cell_lengths_angles_to_cell(cell_lengths, cell_angles)

    def test_read_amber_traj_with_nonorthogonal_cells(self):
        from scipy.io import netcdf_file

        cell_angles = np.array([90.0, 90.0, 120.0])
        with tempfile.TemporaryDirectory() as tmpdir:
            nc_file = os.path.join(tmpdir, "nonorthogonal.nc")
            shutil.copy("amber/02_Heat.nc", nc_file)
            with netcdf_file(nc_file, "a", mmap=False) as f:
                cell_lengths = np.array(f.variables["cell_lengths"][:])
                f.variables["cell_angles"].data[:] = cell_angles

            system = dpdata.LabeledSystem(
                "amber/02_Heat",
                nc_file=nc_file,
                fmt="amber/md",
            )

        cells = system.data["cells"]
        self.assertEqual(system.get_nframes(), cell_lengths.shape[0])
        np.testing.assert_allclose(
            np.linalg.norm(cells, axis=2), cell_lengths, rtol=1e-7, atol=1e-7
        )
        dot_products = np.stack(
            [
                np.sum(cells[:, 1] * cells[:, 2], axis=1),
                np.sum(cells[:, 0] * cells[:, 2], axis=1),
                np.sum(cells[:, 0] * cells[:, 1], axis=1),
            ],
            axis=1,
        )
        computed_angles = np.rad2deg(
            np.arccos(
                dot_products / cell_lengths[:, [1, 0, 0]] / cell_lengths[:, [2, 2, 1]]
            )
        )
        np.testing.assert_allclose(
            computed_angles,
            np.broadcast_to(cell_angles, computed_angles.shape),
            rtol=1e-7,
            atol=1e-7,
        )
        self.assertTrue(np.any(np.abs(cells[:, 1, 0]) > 1e-7))


@unittest.skipIf(
    skip_parmed_related_test, "skip parmed related test. install parmed to fix"
)
class TestAmberMDTarget(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        ll = "amber/corr/low_level"
        ncfile = "amber/corr/rc.nc"
        parmfile = "amber/corr/qmmm.parm7"
        target = ":1"
        self.system_1 = dpdata.LabeledSystem(
            ll,
            nc_file=ncfile,
            parm7_file=parmfile,
            fmt="amber/md",
            use_element_symbols=target,
        )
        self.system_2 = dpdata.LabeledSystem("amber/corr/dp_ll", fmt="deepmd/npy")

        self.places = 5
        self.e_places = 4
        self.f_places = 6
        self.v_places = 6


if __name__ == "__main__":
    unittest.main()
