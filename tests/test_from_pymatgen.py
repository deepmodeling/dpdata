import os
import unittest

from context import dpdata

try:
    from pymatgen.core import Structure  # noqa: F401

    exist_module = True
except Exception:
    exist_module = False


@unittest.skipIf(not exist_module, "skip pymatgen")
class TestFormPytmatgen(unittest.TestCase):
    def setUp(self):
        structure= Structure.from_file(os.path.join('poscars','POSCAR.P42nmc'))
        self.system_1 = dpdata.System(structure,fmt='pymatgen/structure')
        self.system_2 = dpdata.System(os.path.join('poscars','POSCAR.P42nmc'),fmt='poscar')
        self.system_1_cell = self.system_1['cells'][0]
        self.system_2_cell = self.system_2['cells'][0]
        self.system_1_coords = self.system_1['coords'][0]
        self.system_2_coords = self.system_2['coords'][0]
        self.places = 6
    
    def test_all(self):
        self.assertEqual(self.system_1['atom_names'], self.system_2['atom_names'])
        self.assertAlmostEqual(self.system_1['atom_numbs'], self.system_2['atom_numbs'])
        for i in range(len(self.system_1['atom_types'])):
            self.assertEqual(self.system_1['atom_types'][i], self.system_2['atom_types'][i])
        for i in range(self.system_2_cell.shape[0]):
            for j in range(self.system_2_cell.shape[1]):
                self.assertAlmostEqual(self.system_1_cell[i][j], self.system_2_cell[i][j])
        for row in range(self.system_2_coords.shape[0]):
            for col in range(self.system_2_coords.shape[1]):
                self.assertAlmostEqual(self.system_1_coords[row][col], self.system_2_coords[row][col])

if __name__ == "__main__":
    unittest.main()