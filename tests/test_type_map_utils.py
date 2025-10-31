from __future__ import annotations

import unittest

import numpy as np
from context import dpdata

from dpdata.utils import sort_atom_names


class TestSortAtomNames(unittest.TestCase):
    def test_sort_atom_names_type_map(self):
        # Test basic functionality with type_map
        data = {
            "atom_names": ["H", "O"],
            "atom_numbs": [2, 1],
            "atom_types": np.array([1, 0, 0]),
        }
        type_map = ["O", "H"]
        result = sort_atom_names(data, type_map=type_map)
        
        self.assertEqual(result["atom_names"], ["O", "H"])
        self.assertEqual(result["atom_numbs"], [1, 2])
        np.testing.assert_array_equal(result["atom_types"], np.array([0, 1, 1]))
    
    def test_sort_atom_names_type_map_with_zero_atoms(self):
        # Test with type_map that includes elements with zero atoms
        data = {
            "atom_names": ["H", "O"],
            "atom_numbs": [2, 1],
            "atom_types": np.array([1, 0, 0]),
        }
        type_map = ["O", "H", "C"]  # C is not in atom_names but in type_map
        result = sort_atom_names(data, type_map=type_map)
        
        self.assertEqual(result["atom_names"], ["O", "H", "C"])
        self.assertEqual(result["atom_numbs"], [1, 2, 0])
        np.testing.assert_array_equal(result["atom_types"], np.array([0, 1, 1]))
    
    def test_sort_atom_names_type_map_missing_active_types(self):
        # Test that ValueError is raised when active atom types are missing from type_map
        data = {
            "atom_names": ["H", "O"],
            "atom_numbs": [2, 1],  # Both H and O are active (numb > 0)
            "atom_types": np.array([1, 0, 0]),
        }
        type_map = ["H"]  # O is active but missing from type_map
        
        with self.assertRaises(ValueError) as cm:
            sort_atom_names(data, type_map=type_map)
        
        self.assertIn("Active atom types", str(cm.exception))
        self.assertIn("not in provided type_map", str(cm.exception))
        self.assertIn("O", str(cm.exception))
    
    def test_sort_atom_names_without_type_map(self):
        # Test sorting without type_map (alphabetical order)
        data = {
            "atom_names": ["Zn", "O", "H"],
            "atom_numbs": [1, 1, 2],
            "atom_types": np.array([0, 1, 2, 2]),
        }
        result = sort_atom_names(data)
        
        self.assertEqual(result["atom_names"], ["H", "O", "Zn"])
        self.assertEqual(result["atom_numbs"], [2, 1, 1])
        np.testing.assert_array_equal(result["atom_types"], np.array([2, 1, 0, 0]))

    def test_sort_atom_names_with_zero_count_elements_removed(self):
        # Test the case where original elements are A B C, but counts are 0 1 2,
        # which should be able to map to B C (removing A which has count 0)
        data = {
            "atom_names": ["Cl", "O", "C"],
            "atom_numbs": [0, 1, 2],
            "atom_types": np.array([1, 2, 2]),
        }
        type_map = ["O", "C"]  # A is omitted because it has 0 atoms
        result = sort_atom_names(data, type_map=type_map)
        
        self.assertEqual(result["atom_names"], ["O", "C"])
        self.assertEqual(result["atom_numbs"], [1, 2])
        np.testing.assert_array_equal(result["atom_types"], np.array([0, 1, 1]))


if __name__ == "__main__":
    unittest.main()