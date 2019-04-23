import os
import numpy as np
import unittest
from context import dpdata


class TestPWSCFProps :
    def test_atom_names(self) :
        self.assertEqual(self.system.data['atom_names'], ['O','H'])

    def test_atom_numbs(self) :
        self.assertEqual(self.system.data['atom_numbs'], [64,127])

    def test_atom_types(self) :
        for ii in range(0,64) :
            self.assertEqual(self.system.data['atom_types'][ii], 0)
        for ii in range(64,191) :
            self.assertEqual(self.system.data['atom_types'][ii], 1)

    def test_cell(self) :
        ref = 0.52917721067 * 23.5170 * np.eye(3)
        self.assertEqual(self.system.get_nframes(), 2)
        for ff in range(self.system.get_nframes()) :
            for ii in range(3) :
                for jj in range(3) :
                    self.assertEqual(self.system['cells'][ff][ii][jj], ref[ii][jj])

    def test_coord(self) :        
        with open('pwscf.traj/oh-md.pos') as fp :
            lines = fp.read().rstrip('\n').split('\n')
        lines = lines[-191:]
        coords = []
        for ii in lines :
            coords.append([float(jj) for jj in ii.split()])
        bohr2ang = 0.52917721067
        coords = bohr2ang * np.array(coords)
        celll = bohr2ang * 23.5170 
        for ii in range(coords.shape[0]) :
            for jj in range(coords[ii].size) :
                if coords[ii][jj] < 0 :
                    coords[ii][jj] += celll
                elif coords[ii][jj] >= celll :
                    coords[ii][jj] -= celll
                self.assertAlmostEqual(self.system['coords'][-1][ii][jj], coords[ii][jj])


class TestPWSCFTraj(unittest.TestCase, TestPWSCFProps):    

    def setUp(self): 
        self.system = dpdata.System('pwscf.traj/oh-md', fmt = 'pwscf/traj')


class TestPWSCFLabeledTraj(unittest.TestCase, TestPWSCFProps):    

    def setUp(self): 
        self.system = dpdata.LabeledSystem('pwscf.traj/oh-md', fmt = 'pwscf/traj')




if __name__ == '__main__':
    unittest.main()
    
