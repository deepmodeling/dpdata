import unittest
from context import dpdata
from itertools import permutations


class TestTypeMap():
    def test_check_type_map(self):
        # read atom names
        system = dpdata.LabeledSystem(self.fn, fmt=self.fmt)
        atom_names = system['atom_names']
        for type_map in permutations(atom_names, len(atom_names)):
            type_map = list(type_map)
            system.check_type_map(type_map=type_map)
            self.assertEqual(type_map, system['atom_names'])
    
    def test_type_map_is_superset(self):
        system = dpdata.LabeledSystem(self.fn, fmt=self.fmt)
        atom_names = system['atom_names'] + ["X"]
        for type_map in permutations(atom_names, len(atom_names)):
            type_map = list(type_map)
            system = dpdata.LabeledSystem(self.fn, fmt=self.fmt)
            system.check_type_map(type_map=type_map)
            self.assertEqual(type_map, system['atom_names'])


class TestTypeMap1(TestTypeMap, unittest.TestCase):
    def setUp(self):
        self.fn = 'gaussian/methane.gaussianlog'
        self.fmt = 'gaussian/log'


class TestTypeMap2(TestTypeMap, unittest.TestCase):
    def setUp(self):
        self.fn = 'cp2k/cp2k_output'
        self.fmt = 'cp2k/output'


if __name__ == '__main__':
    unittest.main()
