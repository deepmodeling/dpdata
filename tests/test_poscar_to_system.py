import os
import numpy as np
import unittest
from context import dpdata
from poscars.poscar_ref_oh import TestPOSCARoh        

class TestPOSCARCart(unittest.TestCase, TestPOSCARoh):
    
    def setUp(self): 
        self.system = dpdata.System()
        self.system.from_vasp_poscar(os.path.join('poscars', 'POSCAR.oh.c'))

class TestPOSCARDirect(unittest.TestCase, TestPOSCARoh):
    
    def setUp(self): 
        self.system = dpdata.System()
        self.system.from_vasp_poscar(os.path.join('poscars', 'POSCAR.oh.d'))



if __name__ == '__main__':
    unittest.main()
