"""Test that DFTB+ parser handles systems with more than 4 atoms."""

from __future__ import annotations

import os
import sys
import unittest

import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dpdata.formats.dftbplus.output import read_dftb_plus


class TestDftbPlusSixAtoms(unittest.TestCase):
    """Test DFTB+ parser with a 6-atom system (beyond the original 4-atom limit)."""

    def test_six_atoms_parsed_correctly(self):
        """Verify that all 6 atoms are parsed, not truncated to 4."""
        symbols, coord, energy, forces = read_dftb_plus(
            "dftbplus/dftb_pin_6atoms.hsd",
            "dftbplus/detailed_6atoms.out",
        )

        # Should have 6 atoms, not 4
        self.assertEqual(len(symbols), 6)
        self.assertEqual(coord.shape, (6, 3))
        self.assertEqual(forces.shape, (6, 3))

        # Check atom symbols
        expected_symbols = ["C", "H", "O", "N", "S", "Cl"]
        np.testing.assert_array_equal(symbols, expected_symbols)

        # Check coordinates
        expected_coord = np.array(
            [
                [1.0, 0.0, 0.0],
                [1.5, 0.0, 0.0],
                [2.0, 1.0, 0.0],
                [2.5, 0.0, 1.0],
                [3.0, 1.0, 1.0],
                [3.5, 0.5, 0.5],
            ]
        )
        np.testing.assert_array_almost_equal(coord, expected_coord)

        # Check forces
        expected_forces = np.array(
            [
                [0.01, 0.002, 0.003],
                [-0.015, 0.001, 0.002],
                [0.005, -0.008, 0.001],
                [0.003, 0.004, -0.006],
                [-0.002, 0.003, 0.005],
                [0.001, -0.001, -0.002],
            ]
        )
        np.testing.assert_array_almost_equal(forces, expected_forces)

        # Check energy
        self.assertAlmostEqual(energy, -12.3456789012)


if __name__ == "__main__":
    unittest.main()
