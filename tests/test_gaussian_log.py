import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompLabeledSys

class TestGaussianLog :
    def test_atom_names(self) :
        self.assertEqual(self.system.data['atom_names'], ['C','H'])

    def test_atom_numbs(self) :
        self.assertEqual(self.system.data['atom_numbs'], [1, 4])

    def test_atom_types(self) :
        for ii in range(0,1) :
            self.assertEqual(self.system.data['atom_types'][ii], 0)
        for ii in range(1,5) :
            self.assertEqual(self.system.data['atom_types'][ii], 1)

class TestGaussianLoadLog(unittest.TestCase, TestGaussianLog):
    def setUp (self) :
        self.system = dpdata.LabeledSystem('gaussian/methane.gaussianlog', 
                                           fmt = 'gaussian/log')

class TestNonCoveragedGaussianLog :
    def test_atom_names(self) :
        self.assertEqual(self.system.data['atom_names'], [])

    def test_atom_numbs(self) :
        self.assertEqual(self.system.data['atom_numbs'], [])

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
        self.assertEqual(self.system.data['virials'], [])

class TestNonCoveragedGaussianLoadLog(unittest.TestCase, TestNonCoveragedGaussianLog):
    def setUp (self) :
        self.system = dpdata.LabeledSystem('gaussian/noncoveraged.gaussianlog',
                                           fmt = 'gaussian/log')

if __name__ == '__main__':
    unittest.main()
