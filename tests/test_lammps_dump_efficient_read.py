#!/usr/bin/env python3
"""Test efficient frame reading functionality for LAMMPS dump files."""

from __future__ import annotations

import os
import unittest

import numpy as np
from comp_sys import CompSys, IsPBC
from context import dpdata
import dpdata.lammps.dump as dump


class TestLAMMPSDumpEfficientRead(unittest.TestCase, CompSys, IsPBC):
    def setUp(self):
        self.dump_file = os.path.join("poscars", "conf.dump")
        self.type_map = ["O", "H"]
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4
        
        # Set up comparison systems for inherited tests
        # Use the new efficient method as system_1
        self.system_1 = dpdata.System(self.dump_file, fmt="lammps/dump", type_map=self.type_map, f_idx=[0])
        # Use traditional method as system_2
        self.system_2 = dpdata.System(self.dump_file, fmt="lammps/dump", type_map=self.type_map).sub_system([0])

    def test_get_frame_nlines(self):
        """Test frame line count detection."""
        nlines = dump.get_frame_nlines(self.dump_file)
        self.assertEqual(nlines, 11)  # Expected based on file structure

    def test_read_frames_single(self):
        """Test reading a single frame."""
        lines = dump.read_frames(self.dump_file, [1])
        self.assertEqual(len(lines), 11)
        self.assertTrue(lines[0].startswith("ITEM: TIMESTEP"))
        self.assertEqual(lines[1], "1")  # Second frame has timestep 1

    def test_read_frames_multiple(self):
        """Test reading multiple frames."""
        lines = dump.read_frames(self.dump_file, [0, 1])
        self.assertEqual(len(lines), 22)  # 11 lines per frame * 2 frames

    def test_read_frames_out_of_order(self):
        """Test reading frames in non-sequential order."""
        lines1 = dump.read_frames(self.dump_file, [1, 0])
        lines2 = dump.read_frames(self.dump_file, [0, 1])
        self.assertEqual(len(lines1), len(lines2))

    def test_read_frames_empty(self):
        """Test reading with empty frame list."""
        lines = dump.read_frames(self.dump_file, [])
        self.assertEqual(len(lines), 0)

    def test_load_file_with_f_idx(self):
        """Test enhanced load_file with f_idx parameter."""
        # Load specific frame
        lines = dump.load_file(self.dump_file, f_idx=[1])
        self.assertEqual(len(lines), 11)
        
        # Load multiple frames
        lines = dump.load_file(self.dump_file, f_idx=[0, 1])
        self.assertEqual(len(lines), 22)
        
        # Test that f_idx overrides begin/step
        lines = dump.load_file(self.dump_file, begin=1, step=1, f_idx=[0])
        self.assertEqual(len(lines), 11)

    def test_system_with_f_idx(self):
        """Test dpdata.System with f_idx parameter."""
        # Load all frames for comparison
        system_all = dpdata.System(self.dump_file, fmt="lammps/dump", type_map=self.type_map)
        
        # Load only second frame
        system_f1 = dpdata.System(self.dump_file, fmt="lammps/dump", type_map=self.type_map, f_idx=[1])
        
        self.assertEqual(len(system_all.data["coords"]), 2)
        self.assertEqual(len(system_f1.data["coords"]), 1)
        
        # Check that the frame data matches
        np.testing.assert_array_almost_equal(
            system_all.data["coords"][1], 
            system_f1.data["coords"][0]
        )
        np.testing.assert_array_almost_equal(
            system_all.data["cells"][1], 
            system_f1.data["cells"][0]
        )

    def test_load_frames_from_trajectories(self):
        """Test the frames_dict pattern."""
        frames_dict = {
            self.dump_file: [0, 1]
        }
        
        data = dump.load_frames_from_trajectories(frames_dict, type_map=self.type_map)
        
        self.assertIn("coords", data)
        self.assertIn("cells", data)
        self.assertEqual(len(data["coords"]), 2)
        self.assertEqual(len(data["cells"]), 2)

    def test_load_frames_from_trajectories_single(self):
        """Test the frames_dict pattern with single frame."""
        frames_dict = {
            self.dump_file: [1]
        }
        
        data = dump.load_frames_from_trajectories(frames_dict, type_map=self.type_map)
        
        self.assertIn("coords", data)
        self.assertIn("cells", data)
        self.assertEqual(len(data["coords"]), 1)
        self.assertEqual(len(data["cells"]), 1)

    def test_efficiency_comparison(self):
        """Compare efficiency by verifying we get the same results."""
        # Traditional approach: load all then filter
        system_traditional = dpdata.System(self.dump_file, fmt="lammps/dump", type_map=self.type_map)
        filtered_traditional = system_traditional.sub_system([1])
        
        # New efficient approach: load only frame 1
        system_efficient = dpdata.System(self.dump_file, fmt="lammps/dump", type_map=self.type_map, f_idx=[1])
        
        # Results should be identical
        np.testing.assert_array_almost_equal(
            filtered_traditional.data["coords"][0], 
            system_efficient.data["coords"][0]
        )
        np.testing.assert_array_almost_equal(
            filtered_traditional.data["cells"][0], 
            system_efficient.data["cells"][0]
        )

    def setUp_comp_sys(self):
        """Set up comparison systems for inherited tests."""
        pass  # Already set up in setUp


if __name__ == "__main__":
    unittest.main()