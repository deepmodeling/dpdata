import os
import numpy as np
import unittest
from context import dpdata
from poscars.poscar_ref_oh import TestPOSCARoh        

class TestDump(unittest.TestCase, TestPOSCARoh):
    
    def setUp(self): 
        self.system = dpdata.System(os.path.join('poscars', 'conf_unfold.dump'), 
                                    type_map = ['O', 'H'])
        
class TestDump2(unittest.TestCase, TestPOSCARoh):
    
    def setUp(self): 
        self.tmp_system = dpdata.System(os.path.join('poscars', 'conf_unfold.dump'), 
                                        type_map = ['O', 'H'])
        self.system = self.tmp_system.sub_system([1])

    def test_nframes (self) :
        self.assertEqual(self.tmp_system.get_nframes(), 2)
        
        
if __name__ == '__main__':
    unittest.main()
    
