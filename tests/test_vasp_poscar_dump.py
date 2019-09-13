import os
import numpy as np
import unittest
from context import dpdata
from poscars.poscar_ref_oh import TestPOSCARoh        

def myfilecmp(test, f0, f1):
    with open(f0) as fp0 :
        with open(f1) as fp1:
            test.assertTrue(fp0.read() == fp1.read())

class TestPOSCARDump(unittest.TestCase, TestPOSCARoh):
    
    def setUp(self): 
        tmp_system = dpdata.System()
        # tmp_system.from_vasp_poscar(os.path.join('poscars', 'POSCAR.oh.d'))
        tmp_system.from_lammps_lmp(os.path.join('poscars', 'conf.lmp'), type_map = ['O', 'H'])
        tmp_system.to_vasp_poscar('tmp.POSCAR')
        self.system = dpdata.System()
        self.system.from_vasp_poscar('tmp.POSCAR')

class TestPOSCARDump1(unittest.TestCase, TestPOSCARoh):
    
    def setUp(self): 
        tmp_system = dpdata.System()
        tmp_system.from_vasp_poscar(os.path.join('poscars', 'POSCAR.oh.d'))
        # tmp_system.from_lammps_lmp(os.path.join('poscars', 'conf.lmp'), type_map = ['O', 'H'])
        tmp_system.to_vasp_poscar('tmp.POSCAR')
        self.system = dpdata.System()
        self.system.from_vasp_poscar('tmp.POSCAR')

class TestPOSCARSkipZeroAtomNumb(unittest.TestCase) :
    def tearDown(self):
        if os.path.isfile('POSCAR.tmp.1'):
            os.remove('POSCAR.tmp.1')
        if os.path.isfile('POSCAR.tmp.2'):
            os.remove('POSCAR.tmp.2')

    def test_dump_vasp_type_map(self):
        system0 = dpdata.System(os.path.join('poscars', 'POSCAR.oh.d'), fmt = 'vasp/poscar', type_map = ['H', 'O'])
        system0.to_vasp_poscar('POSCAR.tmp.1')
        system1 = dpdata.System(os.path.join('poscars', 'POSCAR.oh.d'), fmt = 'vasp/poscar', type_map = ['C', 'H', 'A', 'O', 'B'])
        system1.to_vasp_poscar('POSCAR.tmp.2')
        myfilecmp(self, 'POSCAR.tmp.1', 'POSCAR.tmp.2')


if __name__ == '__main__':
    unittest.main()
    
