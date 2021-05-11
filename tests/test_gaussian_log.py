import os
import numpy as np
import unittest
from context import dpdata

class TestGaussianLog :
    def test_atom_names(self) :
        self.assertEqual(self.system.data['atom_names'], self.atom_names)

    def test_atom_numbs(self) :
        self.assertEqual(self.system.data['atom_numbs'], self.atom_numbs)
    
    def test_nframes(self):
        self.assertEqual(len(self.system), self.nframes)

    def test_atom_types(self) :
        for ii in range(len(self.atom_types)) :
            self.assertEqual(self.system.data['atom_types'][ii], self.atom_types[ii])

    def test_nopbc(self):
        self.assertEqual(self.system.nopbc, True)

class TestGaussianLoadLog(unittest.TestCase, TestGaussianLog):
    def setUp (self) :
        self.system = dpdata.LabeledSystem('gaussian/methane.gaussianlog', 
                                           fmt = 'gaussian/log')
        self.atom_names = ['C','H']
        self.atom_numbs = [1, 4]
        self.nframes = 1
        self.atom_types = [0, 1, 1, 1, 1]

class TestGaussianLoadLargeForceLog(unittest.TestCase, TestGaussianLog):
    def setUp (self) :
        self.system = dpdata.LabeledSystem('gaussian/largeforce.gaussianlog', 
                                           fmt = 'gaussian/log')
        self.atom_names = ['C','H','O','S']
        self.atom_numbs = [33 , 65, 22, 6]
        self.nframes = 1
        self.atom_types = [0] * 33 + [2] * 22 + [1] * 65 + [3] * 6
    
class TestGaussianLoadMD(unittest.TestCase, TestGaussianLog):
    def setUp (self) :
        self.system = dpdata.LabeledSystem('gaussian/aimd_gaussian_CH4_output', 
                                           fmt = 'gaussian/md')
        self.atom_names = ['C','H']
        self.atom_numbs = [1, 4]
        self.nframes = 22
        self.atom_types = [1, 1, 1, 1, 0]


class TestNonCoveragedGaussianLoadLog(unittest.TestCase, TestGaussianLog):
    def setUp (self) :
        self.system = dpdata.LabeledSystem('gaussian/noncoveraged.gaussianlog',
                                           fmt = 'gaussian/log')
        self.atom_names = []
        self.atom_numbs = []
        self.nframes = 0
    
    def test_atom_types(self) :
        self.assertEqual(self.system.data['atom_types'], [])

    def test_cells(self) :
        self.assertEqual(self.system.data['cells'], [])

    def test_coords(self) :
        self.assertEqual(self.system.data['coords'], [])

    def test_energies(self) :
        self.assertEqual(self.system.data['energies'], [])

    def test_forces(self) :
        self.assertEqual(self.system.data['forces'], [])

    def test_virials(self) :
        self.assertFalse('virials' in self.system.data)

if __name__ == '__main__':
    unittest.main()
