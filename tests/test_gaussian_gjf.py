import os
import unittest
from context import dpdata


class TestDumpGaussianGjf(unittest.TestCase):
    def setUp(self):
        self.system = dpdata.LabeledSystem('gaussian/methane.gaussianlog', 
                                           fmt = 'gaussian/log')
    
    def test_dump_to_gjf(self):
        self.system.to_gaussian_gjf("gaussian/tmp.gjf")
        with open("gaussian/tmp.gjf") as f:
            header = f.readline().strip()
            f.readline()
            title = f.readline().strip()
            f.readline()
            charge, mult = (int(x) for x in f.readline().strip().split())
            atoms = []
            coords = []
            for ii in range(5):
                line = f.readline().strip().split()
                atoms.append(line[0])
                coords.append([float(x) for x in line[1:]])

        self.assertEqual(header, "#force B3LYP/6-31G(d)")
        self.assertEqual(title, self.system.formula)
        self.assertEqual(charge, 0)
        self.assertEqual(mult, 1)
        self.assertEqual(atoms, ['C', 'H', 'H', 'H', 'H'])
        for i in range(self.system['coords'].shape[1]):
            for j in range(3):
                self.assertAlmostEqual(coords[i][j], self.system['coords'][0][i][j])
    
    def tearDown(self):
        if os.path.exists('gaussian/tmp.gjf'):
            os.remove('gaussian/tmp.gjf')