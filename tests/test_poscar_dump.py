import os
import numpy as np
import unittest
from context import system
from poscars.poscar_ref_oh import TestPOSCARoh        

class TestPOSCARDump(unittest.TestCase, TestPOSCARoh):
    
    def setUp(self): 
        tmp_system = system.System()
        # tmp_system.from_vasp_poscar(os.path.join('poscars', 'POSCAR.oh.d'))
        tmp_system.from_lammps_lmp(os.path.join('poscars', 'conf.lmp'), type_map = ['O', 'H'])
        tmp_system.to_vasp_poscar('tmp.POSCAR')
        self.system = system.System()
        self.system.from_vasp_poscar('tmp.POSCAR')

class TestPOSCARDump1(unittest.TestCase, TestPOSCARoh):
    
    def setUp(self): 
        tmp_system = system.System()
        tmp_system.from_vasp_poscar(os.path.join('poscars', 'POSCAR.oh.d'))
        # tmp_system.from_lammps_lmp(os.path.join('poscars', 'conf.lmp'), type_map = ['O', 'H'])
        tmp_system.to_vasp_poscar('tmp.POSCAR')
        self.system = system.System()
        self.system.from_vasp_poscar('tmp.POSCAR')

if __name__ == '__main__':
    unittest.main()
    
