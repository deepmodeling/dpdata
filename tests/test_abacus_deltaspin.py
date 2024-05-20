from __future__ import annotations

import os
import shutil
import unittest

import numpy as np
from context import dpdata

from dpdata.unit import LengthConversion

class TestDeltaSpin():
    """Make a test for the DeltaSpin module in ABACUS
    
    """

    def OUTABACUS_to_system():
        data = dpdata.LabeledSystem('./abacus.scf.deltaspin', fmt = 'abacus/scf')
        print(data)
        
    def system_to_npy():
        data = dpdata.LabeledSystem('./abacus.scf.deltaspin', fmt = 'abacus/scf').to_deepmd_npy('./abacus.scf.deltaspin/')
        print(data)
        



if __name__ == "__main__":
    unittest.main()