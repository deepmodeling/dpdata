import os
import numpy as np
import unittest
from context import dpdata
from poscars.poscar_ref_oh import TestPOSCARoh        

class TestDump(unittest.TestCase, TestPOSCARoh):
    
    def setUp(self): 
        self.system = dpdata.System()
        self.system.from_lammps_dump(os.path.join('poscars', 'conf.dump'), 
                                     type_map = ['O', 'H'])
        
if __name__ == '__main__':
    unittest.main()
    
