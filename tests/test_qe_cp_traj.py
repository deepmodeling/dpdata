import os
import numpy as np
import unittest
from context import dpdata

bohr2ang = dpdata.unit.LengthConversion("bohr", "angstrom").value()

class TestCPTRAJProps :
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
        ref = bohr2ang * 23.5170 * np.eye(3)
        self.assertEqual(self.system.get_nframes(), 2)
        for ff in range(self.system.get_nframes()) :
            for ii in range(3) :
                for jj in range(3) :
                    self.assertEqual(self.system['cells'][ff][ii][jj], ref[ii][jj])

    def test_coord(self) :        
        with open('qe.traj/oh-md.pos') as fp :
            lines = fp.read().rstrip('\n').split('\n')
        lines = lines[-191:]
        coords = []
        for ii in lines :
            coords.append([float(jj) for jj in ii.split()])
        coords = bohr2ang * np.array(coords)
        celll = bohr2ang * 23.5170 
        for ii in range(coords.shape[0]) :
            for jj in range(coords[ii].size) :
                if coords[ii][jj] < 0 :
                    coords[ii][jj] += celll
                elif coords[ii][jj] >= celll :
                    coords[ii][jj] -= celll
                self.assertAlmostEqual(self.system['coords'][-1][ii][jj], coords[ii][jj])


class TestCPTRAJTraj(unittest.TestCase, TestCPTRAJProps):    

    def setUp(self): 
        self.system = dpdata.System('qe.traj/oh-md', fmt = 'qe/cp/traj')


class TestCPTRAJLabeledTraj(unittest.TestCase, TestCPTRAJProps):    

    def setUp(self): 
        self.system = dpdata.LabeledSystem('qe.traj/oh-md', fmt = 'qe/cp/traj')


class TestConverCellDim(unittest.TestCase):    
    def test_case_null(self):
        cell = dpdata.qe.traj.convert_celldm(8, [1, 1, 1])
        ref = np.eye(3)
        for ii in range(3):
            for jj in range(3):
                self.assertAlmostEqual(cell[ii][jj], ref[ii][jj])

class TestCPTRAJTrajCellParamAngstrom(unittest.TestCase):
    def test(self):
        ss = dpdata.System('qe.traj/traj_angs', fmt = 'qe/cp/traj')
        expected_cell = 24 * np.eye(3)
        expected_coord5 = \
            [0.32121855340601E+02,     0.17164728967776E+02,     0.36182628967521E+01,
             0.14321036136090E+01,     0.54322965963498E+01,     0.16038534464162E+02,]
        expected_coord5 = np.array(expected_coord5).reshape([2,3]) * bohr2ang
        self.assertEqual(ss.get_natoms(), 2)
        self.assertEqual(ss.get_nframes(), 6)
        self.assertEqual(ss['atom_names'], ['O','H'])
        self.assertEqual(ss['atom_numbs'], [1,1])
        np.testing.assert_equal(ss['atom_types'], [0,1])
        for ii in range(6):
            np.testing.assert_almost_equal(expected_cell, ss['cells'][ii])
        np.testing.assert_almost_equal(expected_coord5, ss['coords'][5])

class TestCPTRAJTrajCellParamBohr(unittest.TestCase):
    def test(self):
        ss = dpdata.System('qe.traj/traj_bohr', fmt = 'qe/cp/traj')
        expected_cell = 48 * np.eye(3) * bohr2ang
        expected_coord5 = \
            [0.32121855340601E+02,     0.17164728967776E+02,     0.36182628967521E+01,
             0.14321036136090E+01,     0.54322965963498E+01,     0.16038534464162E+02,]
        expected_coord5 = np.array(expected_coord5).reshape([2,3]) * bohr2ang
        self.assertEqual(ss.get_natoms(), 2)
        self.assertEqual(ss.get_nframes(), 6)
        self.assertEqual(ss['atom_names'], ['O','H'])
        self.assertEqual(ss['atom_numbs'], [1,1])
        np.testing.assert_equal(ss['atom_types'], [0,1])
        for ii in range(6):
            np.testing.assert_almost_equal(expected_cell, ss['cells'][ii])
        np.testing.assert_almost_equal(expected_coord5, ss['coords'][5])

class TestCPTRAJTrajCellParamAlat(unittest.TestCase):
    def test(self):
        with self.assertRaises(RuntimeError) as ee:
            ss = dpdata.System('qe.traj/traj_alat', fmt = 'qe/cp/traj')

if __name__ == '__main__':
    unittest.main()
    
