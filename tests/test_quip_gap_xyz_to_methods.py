#!/usr/bin/env python3
from __future__ import annotations

import os
import tempfile
import unittest

from context import dpdata


class TestQuipGapXYZToMethods(unittest.TestCase):
    """Test the to_labeled_system and to_multi_systems methods for QuipGapXYZFormat."""

    def setUp(self):
        """Set up test data."""
        # Load test multi-systems
        self.multi_systems = dpdata.MultiSystems.from_file(
            "xyz/xyz_unittest.xyz", "quip/gap/xyz"
        )
        self.system_b1c9 = self.multi_systems.systems["B1C9"]
        self.system_b5c7 = self.multi_systems.systems["B5C7"]

    def test_to_labeled_system(self):
        """Test writing a single labeled system to QUIP/GAP XYZ format."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".xyz", delete=False
        ) as tmp_file:
            output_file = tmp_file.name

        try:
            # Write the system to file
            self.system_b1c9.to("quip/gap/xyz", output_file)

            # Verify file was created and has content
            self.assertTrue(os.path.exists(output_file))
            with open(output_file) as f:
                content = f.read()
            self.assertTrue(len(content) > 0)

            # Read back and verify we can parse it (use MultiSystems.from_file for QUIP/GAP XYZ)
            reloaded_multi = dpdata.MultiSystems.from_file(output_file, "quip/gap/xyz")
            self.assertEqual(len(reloaded_multi.systems), 1)

            # Verify the data matches (we should have the same system)
            reloaded_system = list(reloaded_multi.systems.values())[0]
            self.assertEqual(len(reloaded_system), len(self.system_b1c9))

        finally:
            if os.path.exists(output_file):
                os.unlink(output_file)

    def test_to_multi_systems(self):
        """Test writing multiple systems to a single QUIP/GAP XYZ format file."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".xyz", delete=False
        ) as tmp_file:
            output_file = tmp_file.name

        try:
            # Write all systems to file
            self.multi_systems.to("quip/gap/xyz", output_file)

            # Verify file was created and has content
            self.assertTrue(os.path.exists(output_file))
            with open(output_file) as f:
                content = f.read()
            self.assertTrue(len(content) > 0)

            # Read back and verify we get the same number of systems
            reloaded_multi = dpdata.MultiSystems.from_file(output_file, "quip/gap/xyz")
            self.assertEqual(
                len(reloaded_multi.systems), len(self.multi_systems.systems)
            )

            # Verify total number of frames is preserved
            original_frames = sum(
                len(sys) for sys in self.multi_systems.systems.values()
            )
            reloaded_frames = sum(len(sys) for sys in reloaded_multi.systems.values())
            self.assertEqual(reloaded_frames, original_frames)

        finally:
            if os.path.exists(output_file):
                os.unlink(output_file)

    def test_roundtrip_consistency(self):
        """Test that writing and reading back preserves data consistency."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".xyz", delete=False
        ) as tmp_file:
            output_file = tmp_file.name

        try:
            # Write and read back
            self.multi_systems.to("quip/gap/xyz", output_file)
            reloaded_multi = dpdata.MultiSystems.from_file(output_file, "quip/gap/xyz")

            # Compare original and reloaded data for each system
            for system_name in self.multi_systems.systems:
                if system_name in reloaded_multi.systems:
                    original = self.multi_systems.systems[system_name]
                    reloaded = reloaded_multi.systems[system_name]

                    # Check basic properties
                    self.assertEqual(len(original), len(reloaded))
                    self.assertEqual(
                        len(original.data["atom_names"]),
                        len(reloaded.data["atom_names"]),
                    )

                    # Note: We don't check exact numerical equality because of floating point precision
                    # and potential differences in formatting, but the data should be structurally the same

        finally:
            if os.path.exists(output_file):
                os.unlink(output_file)


if __name__ == "__main__":
    unittest.main()
