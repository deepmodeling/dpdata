from __future__ import annotations

import os
import shutil
import tempfile
import unittest
import warnings

import numpy as np
from context import dpdata

SOURCE = os.path.join("poscars", "OUTCAR.h2o.md")


class TestOutcarInterrupted(unittest.TestCase):
    """An interrupted VASP run leaves a final ionic step without an energy."""

    def setUp(self):
        with open(SOURCE) as fp:
            self.lines = fp.read().split("\n")
        self.energies = [
            i for i, ll in enumerate(self.lines) if "free  energy   TOTEN" in ll
        ]
        self.tmp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp_dir, ignore_errors=True)

    def _load(self, lines):
        path = os.path.join(self.tmp_dir, "OUTCAR")
        with open(path, "w") as fp:
            fp.write("\n".join(lines))
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            system = dpdata.LabeledSystem(path, fmt="vasp/outcar")
        messages = [
            str(ii.message) for ii in caught if "no energy found" in str(ii.message)
        ]
        return system, messages

    def test_complete_file_is_silent(self):
        system, messages = self._load(self.lines)
        self.assertEqual(system.get_nframes(), 3)
        self.assertEqual(messages, [])

    def test_missing_final_energy_is_reported(self):
        # Cut the file just before the last TOTEN record.
        system, messages = self._load(self.lines[: self.energies[-1]])
        self.assertEqual(system.get_nframes(), 2)
        self.assertEqual(len(messages), 1)
        self.assertIn("frame 3", messages[0])

    def test_frames_before_the_cut_are_intact(self):
        reference = dpdata.LabeledSystem(SOURCE, fmt="vasp/outcar")
        system, _ = self._load(self.lines[: self.energies[-1]])
        np.testing.assert_allclose(
            system["energies"], reference["energies"][:2], atol=1e-10
        )
        np.testing.assert_allclose(
            system["coords"], reference["coords"][:2], atol=1e-10
        )

    def test_energies_stay_numeric(self):
        # The original report was np.savetxt choking on an object array after
        # a None energy was appended.
        system, _ = self._load(self.lines[: self.energies[-1]])
        self.assertEqual(system["energies"].dtype, np.dtype("float64"))


if __name__ == "__main__":
    unittest.main()
