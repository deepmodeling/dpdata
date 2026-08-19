from __future__ import annotations

import io
import os
import unittest

import numpy as np
from comp_sys import CompSys, IsPBC
from context import dpdata

from dpdata.formats.lammps import dump


class CountingStringIO(io.StringIO):
    """Track every supported read strategy and whether it rewinds."""

    def __init__(self, value):
        super().__init__(value)
        self.lines_read = 0
        self.max_position = 0
        self.seek_calls = 0

    def _record_position(self):
        self.max_position = max(self.max_position, self.tell())

    def __next__(self):
        self.lines_read += 1
        value = super().__next__()
        self._record_position()
        return value

    def read(self, *args, **kwargs):
        value = super().read(*args, **kwargs)
        self._record_position()
        return value

    def readline(self, *args, **kwargs):
        value = super().readline(*args, **kwargs)
        self._record_position()
        return value

    def readlines(self, *args, **kwargs):
        value = super().readlines(*args, **kwargs)
        self._record_position()
        return value

    def seek(self, *args, **kwargs):
        self.seek_calls += 1
        return super().seek(*args, **kwargs)


class TestLmpDumpSkip(unittest.TestCase, CompSys, IsPBC):
    def setUp(self):
        self.system_1 = dpdata.System(
            os.path.join("poscars", "conf.5.dump"), type_map=["O", "H"], begin=1, step=2
        )
        self.system_2 = dpdata.System(
            os.path.join("poscars", "conf.5.dump"), type_map=["O", "H"], begin=0, step=1
        ).sub_system(np.arange(1, 5, 2))
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4


class TestLmpDumpFrameSelection(unittest.TestCase):
    def setUp(self):
        self.dump_file = os.path.join("poscars", "conf.5.dump")
        self.type_map = ["O", "H"]

    def test_select_frames_preserves_order_and_duplicates(self):
        all_frames = dpdata.System(
            self.dump_file, fmt="lammps/dump", type_map=self.type_map
        )
        frame_indices = np.array([4, 1, 4])

        selected = dpdata.System(
            self.dump_file,
            fmt="lammps/dump",
            type_map=self.type_map,
            f_idx=frame_indices,
        )
        expected = all_frames.sub_system(frame_indices)

        np.testing.assert_allclose(selected["coords"], expected["coords"])
        np.testing.assert_allclose(selected["cells"], expected["cells"])

    def test_select_single_frame_by_integer(self):
        selected = dpdata.System(
            self.dump_file,
            fmt="lammps/dump",
            type_map=self.type_map,
            f_idx=2,
        )
        expected = dpdata.System(
            self.dump_file, fmt="lammps/dump", type_map=self.type_map
        )[2]

        np.testing.assert_allclose(selected["coords"], expected["coords"])
        np.testing.assert_allclose(selected["cells"], expected["cells"])

    def test_file_object_stops_after_last_target_without_rewinding(self):
        with open(self.dump_file) as fp:
            content = fp.read()
        stream = CountingStringIO(content)

        lines = dump.load_file(stream, f_idx=[1])
        frames = dump.split_traj(lines)

        self.assertEqual(frames[0][1], "1")
        first_23_lines = "".join(content.splitlines(keepends=True)[:23])
        self.assertLessEqual(stream.max_position, len(first_23_lines))
        self.assertLess(stream.max_position, len(content))
        self.assertLessEqual(stream.lines_read, 23)
        self.assertEqual(stream.seek_calls, 0)

    def test_blank_trailer_is_ignored_and_extra_text_is_clamped(self):
        with open(self.dump_file) as fp:
            content = fp.read().rstrip("\n") + "\n\n"

        lines = dump.load_file(io.StringIO(content))
        frames = dump.split_traj([*lines, "# end of run"])
        self.assertEqual(len(frames), 5)
        self.assertFalse(any(not line for frame in frames for line in frame))
        self.assertNotIn("# end of run", frames[-1])

        system = dpdata.System(
            io.StringIO(content), fmt="lammps/dump", type_map=self.type_map
        )
        self.assertEqual(system.get_nframes(), 5)

    def test_invalid_frame_indices(self):
        with self.assertRaisesRegex(ValueError, "must not be empty"):
            dump.load_file(self.dump_file, f_idx=[])
        with self.assertRaisesRegex(ValueError, "non-negative"):
            dump.load_file(self.dump_file, f_idx=[-1])
        with self.assertRaisesRegex(IndexError, "out of range"):
            dump.load_file(self.dump_file, f_idx=[5])
        with self.assertRaisesRegex(TypeError, "only integers"):
            dump.load_file(self.dump_file, f_idx=["a"])
        with self.assertRaisesRegex(TypeError, "only integers"):
            dump.load_file(self.dump_file, f_idx=[True])
        with self.assertRaisesRegex(TypeError, "an iterable of integers"):
            dump.load_file(self.dump_file, f_idx=1.5)

    def test_frame_indices_are_mutually_exclusive_with_slice(self):
        with self.assertRaisesRegex(ValueError, "cannot be combined"):
            dump.load_file(self.dump_file, begin=1, f_idx=[2])
