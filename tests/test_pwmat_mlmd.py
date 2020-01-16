import os
import numpy as np
import unittest
import dpdata


class TestSingleStep(unittest.TestCase):

    def setUp(self):
        self.LabeledSystem1 = dpdata.LabeledSystem(os.path.join('pwmat', 'OUT.MLMD'),\
        fmt='movement' )

    def test_mlmd(self) :

        self.assertEqual(self.LabeledSystem1['energies'], -0.2197270691E+03)
        self.assertEqual(self.LabeledSystem1.get_nframes(), 1)
        self.assertEqual(self.LabeledSystem1.get_natoms(), 5)



if __name__ == '__main__':
    unittest.main()
