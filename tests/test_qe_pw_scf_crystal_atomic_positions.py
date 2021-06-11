import os
import numpy as np
import unittest
from context import dpdata

class TestPWSCFCrystalAtomicPosition:

    def test_coord(self) :
        ref_coord = np.array([[0,0,0], [0, 2.02, 2.02], [2.02, 0, 2.02], [2.02, 2.02, 0]])
        for ii in range(ref_coord.shape[0]) :
            for jj in range(ref_coord.shape[1]) :
                self.assertAlmostEqual(self.system_al.data['coords'][0][ii][jj], ref_coord[ii][jj])

class TestPWSCFLabeledOutput(unittest.TestCase, TestPWSCFCrystalAtomicPosition):

    def setUp(self):
        self.system_al = dpdata.LabeledSystem('qe.scf/Al.out',fmt='qe/pw/scf')

class TestPWSCFLabeledOutputListInput(unittest.TestCase, TestPWSCFCrystalAtomicPosition):

    def setUp(self):
        self.system_al = dpdata.LabeledSystem(['qe.scf/Al.in', 'qe.scf/Al.out'], fmt='qe/pw/scf')

if __name__ == '__main__':
    unittest.main()

