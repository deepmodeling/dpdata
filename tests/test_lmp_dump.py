import os
import numpy as np
import unittest
from context import system
from poscars.poscar_ref_oh import TestPOSCARoh        

class TestLmpDump(unittest.TestCase, TestPOSCARoh):
    
    def setUp(self): 
        tmp_system = system.System()
        tmp_system.from_lammps_lmp(os.path.join('poscars', 'conf.lmp'), 
                                   type_map = ['O', 'H'])
        tmp_system.to_lammps_lmp('tmp.lmp')
        self.system = system.System()
        self.system.from_lammps_lmp('tmp.lmp', 
                                    type_map = ['O', 'H'])


if __name__ == '__main__':
    unittest.main()
    
