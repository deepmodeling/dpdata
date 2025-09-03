#!/usr/bin/env python3
"""Test suite for the string-to-float conversion bug fix in VASP parsing.

This tests the fix for issue #611: "could not convert string to float"
"""

import os
import tempfile
import unittest

from dpdata import System


class TestVaspFloatConversionFix(unittest.TestCase):
    """Test robust error handling for malformed VASP files."""
    
    def setUp(self):
        """Set up temporary files for testing."""
        self.temp_files = []
    
    def tearDown(self):
        """Clean up temporary files."""
        for temp_file in self.temp_files:
            try:
                os.unlink(temp_file)
            except FileNotFoundError:
                pass
    
    def _create_temp_file(self, content, suffix='.poscar'):
        """Create a temporary file with given content."""
        with tempfile.NamedTemporaryFile(mode='w', suffix=suffix, delete=False) as f:
            f.write(content)
            f.flush()
            self.temp_files.append(f.name)
            return f.name
    
    def test_poscar_invalid_scale_factor(self):
        """Test POSCAR parsing with invalid scale factor."""
        poscar_content = """Test
INVALID_SCALE
   1.0000000000000000    0.0000000000000000    0.0000000000000000
   0.0000000000000000    1.0000000000000000    0.0000000000000000
   0.0000000000000000    0.0000000000000000    1.0000000000000000
Mg
1
Cartesian
0.0000000000000000  0.0000000000000000  0.0000000000000000
"""
        temp_file = self._create_temp_file(poscar_content)
        
        with self.assertRaises(ValueError) as cm:
            System(temp_file, fmt='vasp/poscar')
        
        error_msg = str(cm.exception)
        self.assertIn("Failed to parse scale factor from POSCAR line 2", error_msg)
        self.assertIn("INVALID_SCALE", error_msg)
    
    def test_poscar_invalid_cell_vector(self):
        """Test POSCAR parsing with invalid cell vector component."""
        poscar_content = """Test
1.0
   INVALID_CELL    0.0000000000000000    0.0000000000000000
   0.0000000000000000    1.0000000000000000    0.0000000000000000
   0.0000000000000000    0.0000000000000000    1.0000000000000000
Mg
1
Cartesian
0.0000000000000000  0.0000000000000000  0.0000000000000000
"""
        temp_file = self._create_temp_file(poscar_content)
        
        with self.assertRaises(ValueError) as cm:
            System(temp_file, fmt='vasp/poscar')
        
        error_msg = str(cm.exception)
        self.assertIn("Failed to parse cell vectors in POSCAR", error_msg)
        self.assertIn("INVALID_CELL", error_msg)
    
    def test_poscar_invalid_coordinates(self):
        """Test POSCAR parsing with invalid coordinate component."""
        poscar_content = """Test
1.0
   1.0000000000000000    0.0000000000000000    0.0000000000000000
   0.0000000000000000    1.0000000000000000    0.0000000000000000
   0.0000000000000000    0.0000000000000000    1.0000000000000000
Mg
1
Cartesian
INVALID_COORD  0.0000000000000000  0.0000000000000000
"""
        temp_file = self._create_temp_file(poscar_content)
        
        with self.assertRaises(ValueError) as cm:
            System(temp_file, fmt='vasp/poscar')
        
        error_msg = str(cm.exception)
        self.assertIn("Failed to parse coordinates in POSCAR", error_msg)
        self.assertIn("INVALID_COORD", error_msg)
    
    def test_poscar_empty_scale_factor(self):
        """Test POSCAR parsing with empty scale factor."""
        poscar_content = """Test

   1.0000000000000000    0.0000000000000000    0.0000000000000000
   0.0000000000000000    1.0000000000000000    0.0000000000000000
   0.0000000000000000    0.0000000000000000    1.0000000000000000
Mg
1
Cartesian
0.0000000000000000  0.0000000000000000  0.0000000000000000
"""
        temp_file = self._create_temp_file(poscar_content)
        
        with self.assertRaises(ValueError) as cm:
            System(temp_file, fmt='vasp/poscar')
        
        error_msg = str(cm.exception)
        self.assertIn("Failed to parse scale factor from POSCAR line 2", error_msg)
    
    def test_poscar_insufficient_cell_components(self):
        """Test POSCAR parsing with insufficient cell vector components."""
        poscar_content = """Test
1.0
   1.0000000000000000    0.0000000000000000
   0.0000000000000000    1.0000000000000000    0.0000000000000000
   0.0000000000000000    0.0000000000000000    1.0000000000000000
Mg
1
Cartesian
0.0000000000000000  0.0000000000000000  0.0000000000000000
"""
        temp_file = self._create_temp_file(poscar_content)
        
        with self.assertRaises((ValueError, IndexError)):
            System(temp_file, fmt='vasp/poscar')
    
    def test_poscar_valid_file_still_works(self):
        """Test that valid POSCAR files still work correctly."""
        poscar_content = """Test System
1.0
   5.0000000000000000    0.0000000000000000    0.0000000000000000
   0.0000000000000000    5.0000000000000000    0.0000000000000000
   0.0000000000000000    0.0000000000000000    5.0000000000000000
Mg
1
Cartesian
0.0000000000000000  0.0000000000000000  0.0000000000000000
"""
        temp_file = self._create_temp_file(poscar_content)
        
        # This should not raise any exception
        sys = System(temp_file, fmt='vasp/poscar')
        
        # Verify the system was parsed correctly
        self.assertEqual(len(sys), 1)  # 1 frame
        self.assertEqual(sys.get_natoms(), 1)  # 1 atom
        self.assertEqual(sys["atom_names"], ["Mg"])
    
    def test_poscar_special_characters_in_numbers(self):
        """Test POSCAR parsing with special characters that might appear in malformed files."""
        poscar_content = """Test
1.0
   1.0000000000000000    0.0000000000000000    0.0000000000000000
   0.0000000000000000    1.0000000000000000    0.0000000000000000
   0.0000000000000000    0.0000000000000000    1.0000000000000000
Mg
1
Cartesian
0.0000000000000000  0.0000000000000000  *INVALID*
"""
        temp_file = self._create_temp_file(poscar_content)
        
        with self.assertRaises(ValueError) as cm:
            System(temp_file, fmt='vasp/poscar')
        
        error_msg = str(cm.exception)
        self.assertIn("Failed to parse coordinates in POSCAR", error_msg)
        self.assertIn("*INVALID*", error_msg)


if __name__ == '__main__':
    unittest.main()