from __future__ import annotations

import unittest

import numpy as np

from dpdata.amber.md import cell_lengths_angles_to_cell


class TestAmberNonOrthogonalCells(unittest.TestCase):
    def test_orthogonal_cell_conversion(self):
        """Test that orthogonal cells (90° angles) work correctly."""
        # Test case: simple cubic cell with a=10, b=15, c=20, all angles=90°
        cell_lengths = np.array([[10.0, 15.0, 20.0]])
        cell_angles = np.array([[90.0, 90.0, 90.0]])

        expected_cell = np.array(
            [[[10.0, 0.0, 0.0], [0.0, 15.0, 0.0], [0.0, 0.0, 20.0]]]
        )

        result_cell = cell_lengths_angles_to_cell(cell_lengths, cell_angles)

        np.testing.assert_allclose(result_cell, expected_cell, rtol=1e-12, atol=1e-14)

    def test_monoclinic_cell_conversion(self):
        """Test monoclinic cell (beta != 90°, alpha=gamma=90°)."""
        # Test case: monoclinic cell with a=10, b=15, c=20, alpha=90°, beta=120°, gamma=90°
        cell_lengths = np.array([[10.0, 15.0, 20.0]])
        cell_angles = np.array([[90.0, 120.0, 90.0]])

        # Expected result:
        # v1 = [10, 0, 0]
        # v2 = [0, 15, 0] (gamma=90°)
        # v3 = [20*cos(120°), 0, 20*sin(120°)] = [-10, 0, 17.32...]
        cos_120 = np.cos(np.deg2rad(120.0))  # -0.5
        sin_120 = np.sin(np.deg2rad(120.0))  # sqrt(3)/2

        expected_cell = np.array(
            [
                [
                    [10.0, 0.0, 0.0],
                    [0.0, 15.0, 0.0],
                    [20.0 * cos_120, 0.0, 20.0 * sin_120],
                ]
            ]
        )

        result_cell = cell_lengths_angles_to_cell(cell_lengths, cell_angles)

        np.testing.assert_allclose(result_cell, expected_cell, rtol=1e-12, atol=1e-14)

    def test_hexagonal_cell_conversion(self):
        """Test hexagonal cell (gamma=120°, alpha=beta=90°)."""
        # Test case: hexagonal cell with a=10, b=10, c=15, alpha=90°, beta=90°, gamma=120°
        cell_lengths = np.array([[10.0, 10.0, 15.0]])
        cell_angles = np.array([[90.0, 90.0, 120.0]])

        # Expected result:
        # v1 = [10, 0, 0]
        # v2 = [10*cos(120°), 10*sin(120°), 0] = [-5, 8.66..., 0]
        # v3 = [0, 0, 15] (alpha=beta=90°)
        cos_120 = np.cos(np.deg2rad(120.0))  # -0.5
        sin_120 = np.sin(np.deg2rad(120.0))  # sqrt(3)/2

        expected_cell = np.array(
            [
                [
                    [10.0, 0.0, 0.0],
                    [10.0 * cos_120, 10.0 * sin_120, 0.0],
                    [0.0, 0.0, 15.0],
                ]
            ]
        )

        result_cell = cell_lengths_angles_to_cell(cell_lengths, cell_angles)

        np.testing.assert_allclose(result_cell, expected_cell, rtol=1e-12, atol=1e-14)

    def test_triclinic_cell_conversion(self):
        """Test triclinic cell (all angles != 90°)."""
        # Test case: triclinic cell with a=8, b=10, c=12, alpha=70°, beta=80°, gamma=110°
        cell_lengths = np.array([[8.0, 10.0, 12.0]])
        cell_angles = np.array([[70.0, 80.0, 110.0]])

        result_cell = cell_lengths_angles_to_cell(cell_lengths, cell_angles)

        # Check that the result has the right shape
        self.assertEqual(result_cell.shape, (1, 3, 3))

        # Check that the cell vectors have the correct lengths
        computed_lengths = np.linalg.norm(result_cell[0], axis=1)
        expected_lengths = np.array([8.0, 10.0, 12.0])
        np.testing.assert_allclose(computed_lengths, expected_lengths, rtol=1e-12)

        # Check that the angles between vectors are correct
        v1, v2, v3 = result_cell[0]

        # Angle between v2 and v3 should be alpha (70°)
        cos_alpha = np.dot(v2, v3) / (np.linalg.norm(v2) * np.linalg.norm(v3))
        alpha_computed = np.rad2deg(np.arccos(cos_alpha))
        np.testing.assert_allclose(alpha_computed, 70.0, rtol=1e-8)

        # Angle between v1 and v3 should be beta (80°)
        cos_beta = np.dot(v1, v3) / (np.linalg.norm(v1) * np.linalg.norm(v3))
        beta_computed = np.rad2deg(np.arccos(cos_beta))
        np.testing.assert_allclose(beta_computed, 80.0, rtol=1e-8)

        # Angle between v1 and v2 should be gamma (110°)
        cos_gamma = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        gamma_computed = np.rad2deg(np.arccos(cos_gamma))
        np.testing.assert_allclose(gamma_computed, 110.0, rtol=1e-8)

    def test_extreme_angles_case(self):
        """Test edge case with angles very far from 90°."""
        cell_lengths = np.array([[5.0, 8.0, 12.0]])
        cell_angles = np.array([[60.0, 70.0, 130.0]])  # all far from 90°

        # Should work without error
        result = cell_lengths_angles_to_cell(cell_lengths, cell_angles)
        self.assertEqual(result.shape, (1, 3, 3))

        # Verify the lengths are preserved
        computed_lengths = np.linalg.norm(result[0], axis=1)
        expected_lengths = np.array([5.0, 8.0, 12.0])
        np.testing.assert_allclose(computed_lengths, expected_lengths, rtol=1e-10)

    def test_multiple_frames(self):
        """Test that multiple frames are handled correctly."""
        # Test case: 3 frames with different cell parameters
        cell_lengths = np.array(
            [
                [10.0, 10.0, 10.0],  # cubic
                [8.0, 12.0, 15.0],  # orthorhombic
                [10.0, 10.0, 12.0],
            ]
        )  # hexagonal-like
        cell_angles = np.array(
            [[90.0, 90.0, 90.0], [90.0, 90.0, 90.0], [90.0, 90.0, 120.0]]
        )

        result_cell = cell_lengths_angles_to_cell(cell_lengths, cell_angles)

        # Check shape
        self.assertEqual(result_cell.shape, (3, 3, 3))

        # Check first frame (cubic)
        expected_frame1 = np.array(
            [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
        )
        np.testing.assert_allclose(
            result_cell[0], expected_frame1, rtol=1e-12, atol=1e-14
        )

        # Check third frame (hexagonal-like)
        cos_120 = np.cos(np.deg2rad(120.0))  # -0.5
        sin_120 = np.sin(np.deg2rad(120.0))  # sqrt(3)/2
        expected_frame3 = np.array(
            [[10.0, 0.0, 0.0], [10.0 * cos_120, 10.0 * sin_120, 0.0], [0.0, 0.0, 12.0]]
        )
        np.testing.assert_allclose(
            result_cell[2], expected_frame3, rtol=1e-12, atol=1e-14
        )


if __name__ == "__main__":
    unittest.main()
