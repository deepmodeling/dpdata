import os
import numpy as np
import unittest
from context import dpdata
from poscars.poscar_ref_oh import TestPOSCARoh        

class TestLmpDump(unittest.TestCase, TestPOSCARoh):
    
    def setUp(self): 
        tmp_system = dpdata.System(os.path.join('poscars', 'conf.lmp'), 
                                   type_map = ['O', 'H'])
        tmp_system.to_lammps_lmp('tmp.lmp')
        self.system = dpdata.System()
        self.system.from_lammps_lmp('tmp.lmp', 
                                    type_map = ['O', 'H'])

class TestToFunc(unittest.TestCase, TestPOSCARoh):
    
    def setUp(self): 
        tmp_system = dpdata.System(os.path.join('poscars', 'conf.lmp'), 
                                   type_map = ['O', 'H'])
        tmp_system.to('lammps/lmp', 'tmp.lmp')
        self.system = dpdata.System()
        self.system.from_fmt('tmp.lmp', fmt='lammps/lmp',
                             type_map = ['O', 'H'])

if __name__ == '__main__':
    unittest.main()
    
