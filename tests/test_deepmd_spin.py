from __future__ import annotations

import os
import shutil
import unittest

from context import dpdata


class TestDeepmdReadSpinNPY(unittest.TestCase):
    def setUp(self):
        self.tmp_save_path = "tmp.deepmd.spin/dump-tmp"

    def tearDown(self):
        if os.path.exists(self.tmp_save_path):
            shutil.rmtree(self.tmp_save_path)

    def check_Fe16(self, system):
        self.assertTrue("spins" in system.data)
        self.assertTrue("force_mags" in system.data)
        self.assertEqual(system.data["spins"].shape, (2, 16, 3))
        self.assertEqual(system.data["force_mags"].shape, (2, 16, 3))

    def test_read_spin_npy(self):
        system = dpdata.LabeledSystem("tmp.deepmd.spin/Fe16-npy", fmt="deepmd/npy")
        self.check_Fe16(system)

        system.to("deepmd/npy", self.tmp_save_path)
        self.assertTrue(
            os.path.isfile(os.path.join(self.tmp_save_path, "set.000/spin.npy"))
        )
        self.assertTrue(
            os.path.isfile(os.path.join(self.tmp_save_path, "set.000/force_mag.npy"))
        )

    def test_read_spin_raw(self):
        system = dpdata.LabeledSystem("tmp.deepmd.spin/Fe16-raw", fmt="deepmd/raw")
        self.check_Fe16(system)

        system.to("deepmd/raw", self.tmp_save_path)
        self.assertTrue(os.path.isfile(os.path.join(self.tmp_save_path, "spin.raw")))
        self.assertTrue(
            os.path.isfile(os.path.join(self.tmp_save_path, "force_mag.raw"))
        )


if __name__ == "__main__":
    unittest.main()
