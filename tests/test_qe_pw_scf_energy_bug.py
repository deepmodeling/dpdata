import os
import numpy as np
import unittest
from context import dpdata

class TestPWSCFSinglePointEnergy:

    def test_energy(self) :
        ref_energy = -296.08379065679094669
        self.assertAlmostEqual(self.system_al.data['energies'][0], ref_energy)

class TestPWSCFLabeledOutput(unittest.TestCase, TestPWSCFSinglePointEnergy):

    def setUp(self):
        self.system_al = dpdata.LabeledSystem('qe.scf/Al.out',fmt='qe/pw/scf')

class TestPWSCFLabeledOutputListInput(unittest.TestCase, TestPWSCFSinglePointEnergy):

    def setUp(self):
        self.system_al = dpdata.LabeledSystem(['qe.scf/Al.in', 'qe.scf/Al.out'], fmt='qe/pw/scf')

if __name__ == '__main__':
    unittest.main()

