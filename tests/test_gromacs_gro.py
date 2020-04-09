import os
import numpy as np
import unittest
from context import dpdata

class TestGromacsGro(unittest.TestCase):
    def test_read_file(self):
        system = dpdata.System('gromacs/1h.gro')
        self.assertEqual(system['atom_names'], ['H', 'O'])
        self.assertEqual(system['atom_numbs'], [6, 3])
        for cc,ii in enumerate([1, 0, 0, 1, 0, 0, 1, 0, 0]):
            self.assertEqual(system['atom_types'][cc], ii)
        self.assertEqual(len(system['cells']), 1)
        self.assertEqual(len(system['coords']), 1)
        for ii in range(3):
            for jj in range(3):
                if ii != jj:
                    self.assertAlmostEqual(system['cells'][0][ii][jj], 0)
        self.assertAlmostEqual(system['cells'][0][0][0], 7.822838765564372)
        self.assertAlmostEqual(system['cells'][0][1][1], 7.353572647182051)
        self.assertAlmostEqual(system['cells'][0][2][2], 9.036518515423753)
        self.assertAlmostEqual(system['coords'][0][8][0], 7.43)
        self.assertAlmostEqual(system['coords'][0][8][1], 5.12)
        self.assertAlmostEqual(system['coords'][0][8][2], 3.36)

    def test_read_file_tri(self):
        system = dpdata.System('gromacs/1h.tri.gro')
        self.assertEqual(system['atom_names'], ['H', 'O'])
        self.assertEqual(system['atom_numbs'], [6, 3])
        for cc,ii in enumerate([1, 0, 0, 1, 0, 0, 1, 0, 0]):
            self.assertEqual(system['atom_types'][cc], ii)
        self.assertEqual(len(system['cells']), 1)
        self.assertEqual(len(system['coords']), 1)
        count = 0
        for ii in range(3):
            for jj in range(3):
                if ii != jj:
                    self.assertAlmostEqual(system['cells'][0][ii][jj], count)
                    count += 1
        self.assertAlmostEqual(system['cells'][0][0][0], 7.822838765564372)
        self.assertAlmostEqual(system['cells'][0][1][1], 7.353572647182051)
        self.assertAlmostEqual(system['cells'][0][2][2], 9.036518515423753)
        self.assertAlmostEqual(system['coords'][0][8][0], 7.43)
        self.assertAlmostEqual(system['coords'][0][8][1], 5.12)
        self.assertAlmostEqual(system['coords'][0][8][2], 3.36)
        system.to('vasp/poscar', 'POSCAR')
