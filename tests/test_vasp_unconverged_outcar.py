import os
import numpy as np
import unittest
from context import dpdata


class TestSingleStep(unittest.TestCase):

    def setUp(self):
        self.LabeledSystem1 = dpdata.LabeledSystem(os.path.join('poscars', 'OUTCAR.ch4.unconverged'),\
        fmt='outcar' )

        self.LabeledSystem2 = dpdata.LabeledSystem(os.path.join('poscars', 'OUTCAR.ch4.1step'),\
        fmt='outcar' )

    def test_unconverged(self) :

        self.assertEqual(self.LabeledSystem1['energies'], -23.94708651)
        self.assertEqual(self.LabeledSystem1.get_nframes(), 1)
        self.assertEqual(self.LabeledSystem1.get_natoms(), 5)
    def test_single_step(self) :
        self.assertEqual(self.LabeledSystem2.get_nframes(), 0)



if __name__ == '__main__':
    unittest.main()