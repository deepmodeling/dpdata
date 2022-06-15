import os
import numpy as np
import unittest
from context import dpdata
from test_vasp_poscar_dump import myfilecmp


class TestStruDump(unittest.TestCase):
    def setUp(self):
        self.system_ch4 = dpdata.System("abacus.scf/STRU.ch4", fmt="stru")

    def test_dump_stru(self):
        self.system_ch4.to("stru", "STRU_tmp", mass = [12, 1], pp_file = ["C.upf", "H.upf"], numerical_orbital = ["C.orb", "H.orb"], numerical_descriptor = "jle.orb")
        myfilecmp(self, "abacus.scf/stru_test", "STRU_tmp")
    
    def tearDown(self):
        if os.path.isfile('STRU_tmp'):
            os.remove('STRU_tmp')

if __name__ == '__main__':
    unittest.main()
    