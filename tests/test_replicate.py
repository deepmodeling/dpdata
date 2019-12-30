import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompSys, IsPBC

class TestReplicate123(unittest.TestCase, CompSys, IsPBC):
    def setUp (self) :
        system_1_origin = dpdata.System('poscars/POSCAR.SiC',fmt='vasp/poscar')
        self.system_1 = system_1_origin.replicate((1,2,3,))
        self.system_2 = dpdata.System('poscars/POSCAR.SiC.replicate123',fmt='vasp/poscar')
        self.places = 6

class TestReplicate123_not_change_origin(unittest.TestCase, CompSys, IsPBC):
    def setUp (self) :
        self.system_1 = dpdata.System('poscars/POSCAR.SiC',fmt='vasp/poscar')
        self.system_1.replicate((1,2,3,))
        self.system_2 = dpdata.System('poscars/POSCAR.SiC',fmt='vasp/poscar')
        self.places = 6

if __name__ == '__main__':
    unittest.main()
