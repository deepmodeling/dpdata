from __future__ import annotations

import os
import shutil
import tempfile
import unittest

import numpy as np
from context import dpdata

HYBRID = os.path.join("poscars", "OUTCAR.ch4.ml")


class TestVaspOutcarMLOnly(unittest.TestCase):
    """An ``ML_ISTART = 2`` run writes no ab initio energy at all.

    ``OUTCAR.ch4.ml`` is a hybrid run, so it still carries
    ``free  energy   TOTEN`` records. Dropping them yields the pure-inference
    layout of #828, where the ab initio token never appears.
    """

    def setUp(self):
        with open(HYBRID) as fp:
            lines = fp.read().split("\n")
        self.pure_ml = [ii for ii in lines if "free  energy   TOTEN" not in ii]
        self.tmp_dir = tempfile.mkdtemp()
        self.path = os.path.join(self.tmp_dir, "OUTCAR")
        with open(self.path, "w") as fp:
            fp.write("\n".join(self.pure_ml))

    def tearDown(self):
        shutil.rmtree(self.tmp_dir, ignore_errors=True)

    def test_all_ml_frames_are_read(self):
        # Before the fix the first block ran to the end of file, so
        # analyze_block returned the first ML energy and the loop stopped:
        # one frame instead of ten.
        system = dpdata.LabeledSystem(self.path, fmt="vasp/outcar", ml=True)
        self.assertEqual(system.get_nframes(), 10)

    def test_matches_the_hybrid_run(self):
        system = dpdata.LabeledSystem(self.path, fmt="vasp/outcar", ml=True)
        reference = dpdata.LabeledSystem(HYBRID, fmt="vasp/outcar", ml=True)
        self.assertEqual(system.get_nframes(), reference.get_nframes())
        np.testing.assert_allclose(system["energies"], reference["energies"])
        np.testing.assert_allclose(system["coords"], reference["coords"])
        np.testing.assert_allclose(system["forces"], reference["forces"])
        self.assertEqual(system["atom_names"], reference["atom_names"])

    def test_system_info_still_parsed(self):
        # The header must still be inside the first block once that block
        # ends at the first ML energy instead of the end of file.
        system = dpdata.LabeledSystem(self.path, fmt="vasp/outcar", ml=True)
        self.assertEqual(system["atom_names"], ["H", "C"])
        self.assertEqual(list(system["atom_types"]), [0, 0, 0, 0, 1])

    def test_hybrid_run_is_unaffected(self):
        self.assertEqual(
            dpdata.LabeledSystem(HYBRID, fmt="vasp/outcar", ml=True).get_nframes(), 10
        )
        self.assertEqual(
            dpdata.LabeledSystem(HYBRID, fmt="vasp/outcar", ml=False).get_nframes(), 4
        )


if __name__ == "__main__":
    unittest.main()
