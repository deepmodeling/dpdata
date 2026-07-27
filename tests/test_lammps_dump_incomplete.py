from __future__ import annotations

import os
import shutil
import tempfile
import unittest
import warnings

import numpy as np
from context import dpdata

SOURCE = os.path.join("poscars", "conf.5.dump")


class TestLmpDumpIncomplete(unittest.TestCase):
    """A dump file may be truncated mid-frame or carry stray comment lines."""

    def setUp(self):
        with open(SOURCE) as fp:
            self.lines = fp.read().split("\n")
        self.tmp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp_dir, ignore_errors=True)

    def _write(self, name, lines):
        path = os.path.join(self.tmp_dir, name)
        with open(path, "w") as fp:
            fp.write("\n".join(lines))
        return path

    def _load(self, path, **kwargs):
        return dpdata.System(path, fmt="lammps/dump", type_map=["O", "H"], **kwargs)

    def test_complete_file_is_unchanged(self):
        self.assertEqual(self._load(SOURCE).get_nframes(), 5)

    def test_truncated_last_frame_is_skipped(self):
        # A run killed while writing leaves the final frame without its atom
        # lines; the earlier frames are still usable.
        path = self._write("truncated.dump", self.lines[:-3])
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            system = self._load(path)
        self.assertEqual(system.get_nframes(), 4)
        self.assertTrue(
            any("incomplete frame 4" in str(ii.message) for ii in caught),
            "dropping a frame must be reported",
        )

    def test_truncated_frame_does_not_shift_later_frames(self):
        # Frames are sliced at their own markers, so a short frame in the
        # middle must not push the remaining frames out of alignment.
        starts = [i for i, ll in enumerate(self.lines) if "ITEM: TIMESTEP" in ll]
        damaged = self.lines[: starts[1] + 4] + self.lines[starts[2] :]
        path = self._write("short_middle.dump", damaged)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            system = self._load(path)
        reference = self._load(SOURCE)
        self.assertEqual(system.get_nframes(), 4)
        self.assertTrue(any("incomplete frame 1" in str(ii.message) for ii in caught))
        # Frame 1 is the damaged one; frames 2..4 must survive intact.
        np.testing.assert_allclose(
            system["coords"][1:], reference["coords"][2:], atol=1e-10
        )

    def test_comment_lines_are_ignored(self):
        atoms = [i for i, ll in enumerate(self.lines) if ll.startswith("ITEM: ATOMS")]
        annotated = list(self.lines)
        for offset, idx in enumerate(atoms):
            annotated.insert(idx + 1 + offset, "##### restart marker")
        path = self._write("commented.dump", annotated)
        system = self._load(path)
        reference = self._load(SOURCE)
        self.assertEqual(system.get_nframes(), reference.get_nframes())
        np.testing.assert_allclose(system["coords"], reference["coords"], atol=1e-10)

    def test_no_usable_frame_raises(self):
        starts = [i for i, ll in enumerate(self.lines) if "ITEM: TIMESTEP" in ll]
        path = self._write("headers_only.dump", self.lines[: starts[0] + 4])
        with self.assertRaisesRegex(RuntimeError, "no complete frame"):
            self._load(path)

    def test_not_a_dump_file_raises(self):
        path = self._write("garbage.dump", ["not a dump file", "1 2 3"])
        with self.assertRaisesRegex(RuntimeError, "not a LAMMPS dump file"):
            self._load(path)


if __name__ == "__main__":
    unittest.main()
