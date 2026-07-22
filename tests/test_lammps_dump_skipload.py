from __future__ import annotations

import io
import os
import unittest

import numpy as np
from comp_sys import CompSys, IsPBC
from context import dpdata

from dpdata.formats.lammps import dump


class CountingStringIO(io.StringIO):
    """Track how many lines the trajectory reader consumes."""

    def __init__(self, value):
        super().__init__(value)
        self.lines_read = 0

    def __next__(self):
        self.lines_read += 1
        return super().__next__()


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

    def test_file_object_is_read_once_and_stops_after_last_target(self):
        with open(self.dump_file) as fp:
            content = fp.read()
        stream = CountingStringIO(content)

        lines = dump.load_file(stream, f_idx=[1])
        frames = dump.split_traj(lines)

        self.assertEqual(frames[0][1], "1")
        self.assertLess(stream.lines_read, len(content.splitlines()))

    def test_invalid_frame_indices(self):
        with self.assertRaisesRegex(ValueError, "must not be empty"):
            dump.load_file(self.dump_file, f_idx=[])
        with self.assertRaisesRegex(ValueError, "non-negative"):
            dump.load_file(self.dump_file, f_idx=[-1])
        with self.assertRaisesRegex(IndexError, "out of range"):
            dump.load_file(self.dump_file, f_idx=[5])

    def test_frame_indices_are_mutually_exclusive_with_slice(self):
        with self.assertRaisesRegex(ValueError, "cannot be combined"):
            dump.load_file(self.dump_file, begin=1, f_idx=[2])
