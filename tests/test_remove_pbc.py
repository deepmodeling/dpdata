import os
import numpy as np
import unittest
from context import dpdata

class TestRemovePBC(unittest.TestCase):

    def test_remove(self):
        coords = np.array([[[-1, -1, 2], [-1,-1,-3], [-1,-1, 7]], 
                           [[ 3, -1, 3], [-1,-1, 3], [ 7,-1, 3]]], dtype = float)
        cogs = np.average(coords, axis = 1)        
        data = {'atom_names' : ['A', 'B'],
                'atom_numbs' : [1, 2],
                'atom_types' : np.array([1, 0, 1], dtype = int),
                'orig': np.array([0, 0, 0]),
                'coords': coords,
                'cells': np.random.random([2, 3, 3]),
        }
        sys = dpdata.System(data = data)
        proct = 9.0
        
        mol_size = np.array([5, 4], dtype = float)
        cell_size = (mol_size + proct) * 2.0

        sys.remove_pbc(proct)

        for ff in range(2):
            ref = cell_size[ff] * np.eye(3)
            for ii in range(3):
                for jj in range(3):
                    self.assertAlmostEqual(sys['cells'][ff][ii][jj], ref[ii][jj], msg = '%d %d %d' %(ff, ii, jj))
            dists = []
            for ii in range(sys.get_natoms()):
                for jj in range(3):
                    dists.append(np.abs(sys['coords'][ff][ii][jj]))
                    dists.append(np.abs(sys['cells'][ff][jj][jj] - sys['coords'][ff][ii][jj]))
            self.assertAlmostEqual(np.min(dists), proct)
