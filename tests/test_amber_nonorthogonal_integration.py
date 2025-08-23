from __future__ import annotations

import unittest

import numpy as np

from dpdata.amber.md import cell_lengths_angles_to_cell


class TestAmberNonOrthogonalIntegration(unittest.TestCase):
    def test_no_runtime_error_for_nonorthogonal_cells(self):
        """Test that non-orthogonal cells no longer raise RuntimeError."""
        # Test case that would have previously raised RuntimeError
        cell_lengths = np.array([[10.0, 10.0, 15.0]])
        cell_angles = np.array([[90.0, 90.0, 120.0]])  # gamma != 90
        
        # This should NOT raise a RuntimeError anymore
        try:
            result = cell_lengths_angles_to_cell(cell_lengths, cell_angles)
            # Should return a valid 3x3 cell matrix
            self.assertEqual(result.shape, (1, 3, 3))
        except RuntimeError as e:
            if "Unsupported cells" in str(e):
                self.fail("cell_lengths_angles_to_cell should support non-orthogonal cells")
            else:
                raise  # Re-raise if it's a different RuntimeError

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

    def test_backwards_compatibility_orthogonal(self):
        """Test that orthogonal cells still produce the same result as before."""
        cell_lengths = np.array([[10.0, 15.0, 20.0]])
        cell_angles = np.array([[90.0, 90.0, 90.0]])
        
        result = cell_lengths_angles_to_cell(cell_lengths, cell_angles)
        
        # Should produce the same diagonal matrix as the old implementation
        expected = np.array([[[10.0, 0.0, 0.0],
                             [0.0, 15.0, 0.0],
                             [0.0, 0.0, 20.0]]])
        
        np.testing.assert_allclose(result, expected, rtol=1e-12, atol=1e-14)


if __name__ == "__main__":
    unittest.main()