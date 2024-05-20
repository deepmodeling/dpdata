from __future__ import annotations

import unittest

from context import dpdata

class TestDeltaSpin(unittest.TestCase):
    """Make a test for the DeltaSpin module in ABACUS
    
    """

    def OUTABACUS_to_system(self):
        data = dpdata.LabeledSystem('./abacus.scf.deltaspin', fmt = 'abacus/scf')
        print(data)
        
    def system_to_npy(self):
        data = dpdata.LabeledSystem('./abacus.scf.deltaspin', fmt = 'abacus/scf').to_deepmd_npy('./abacus.scf.deltaspin/')
        print(data)


if __name__ == "__main__":
    unittest.main()
