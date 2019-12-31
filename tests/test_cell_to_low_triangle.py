import os
import numpy as np
import unittest
from context import dpdata

class TestCellToLowTriangle(unittest.TestCase):
    def test_func1(self):
        cell_1 = dpdata.cp2k.cell.cell_to_low_triangle(6,6,6,np.pi*1/2, np.pi*1/2, np.pi*1/2)
        cell_2 = np.asarray([[6,0,0],[0,6,0],[0,0,6]])
        for ii in range(3):
            for jj in range(3):
                self.assertAlmostEqual(cell_1[ii,jj], cell_2[ii,jj], places=6)

    def test_func2(self):
        cell_1 = dpdata.cp2k.cell.cell_to_low_triangle(6,6,6,np.pi*1/3, np.pi*1/3, np.pi*1/3)
        cell_2 = np.asarray([
            [6,0,0],
            [3,3*np.sqrt(3),0],
            [3,np.sqrt(3),2*np.sqrt(6)]])
        for ii in range(3):
            for jj in range(3):
                self.assertAlmostEqual(cell_1[ii,jj], cell_2[ii,jj], places=6)

    def test_func3(self):
        cell_1 = dpdata.cp2k.cell.cell_to_low_triangle(6,7,8,np.pi*133/180,np.pi*84/180,np.pi*69/180)
        cell_2 = np.asarray([[ 6.0, 0.0, 0.0],
            [ 2.5085757,  6.535063 ,  0.0],
            [ 0.8362277, -6.1651506,  5.0290794]], dtype='float32')
        for ii in range(3):
            for jj in range(3):
                self.assertAlmostEqual(cell_1[ii,jj], cell_2[ii,jj], places=6)

    def test_func4(self):
        with self.assertRaises(Exception) as c:
            dpdata.cp2k.cell.cell_to_low_triangle(0.1,6,6,np.pi*1/2,np.pi*1/2,np.pi*1/2)
        self.assertTrue("A==0.1" in str(c.exception))

    def test_func5(self):
        with self.assertRaises(Exception) as c:
            dpdata.cp2k.cell.cell_to_low_triangle(6,6,6,np.pi*3/180,np.pi*1/2,np.pi*1/2)
        self.assertTrue("alpha" in str(c.exception))

    def test_func6(self):
        with self.assertRaises(Exception) as c:
            dpdata.cp2k.cell.cell_to_low_triangle(6,7,8,np.pi*153/180,np.pi*84/180,np.pi*69/180)
        self.assertTrue("lz^2" in str(c.exception))

if __name__ == '__main__':
    unittest.main()
