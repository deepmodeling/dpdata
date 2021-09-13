import os
import numpy as np
import unittest
from context import dpdata

class TestPOSCARCart(unittest.TestCase):
    
    def setUp(self): 
        self.system = dpdata.System()
        self.system.from_pymatgen_molecule(os.path.join('pymatgen', 'FA-001.xyz'))

    def test_to_molecule(self):
        mols = self.system.to_pymatgen_molecule()
        print(mols[-1])



if __name__ == '__main__':
    unittest.main()

