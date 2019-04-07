import os
import numpy as np
import unittest
from context import system
from poscars.poscar_ref_oh import TestPOSCARoh        

class TestPOSCARCart(unittest.TestCase, TestPOSCARoh):
    
    def setUp(self): 
        self.system = system.System()
        self.system.from_vasp_poscar(os.path.join('poscars', 'POSCAR.oh.c'))

class TestPOSCARDirect(unittest.TestCase, TestPOSCARoh):
    
    def setUp(self): 
        self.system = system.System()
        self.system.from_vasp_poscar(os.path.join('poscars', 'POSCAR.oh.d'))



if __name__ == '__main__':
    unittest.main()
