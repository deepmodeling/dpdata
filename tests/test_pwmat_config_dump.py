import os
import numpy as np
import unittest
import dpdata
from pwmat.config_ref_oh import Testconfigoh       

def myfilecmp(test, f0, f1):
    with open(f0) as fp0 :
        with open(f1) as fp1:
            test.assertTrue(fp0.read() == fp1.read())

class TestatomconfigDump(unittest.TestCase, Testconfigoh):
    
    def setUp(self):
        tmp_system = dpdata.System()
        tmp_system.from_lammps_lmp(os.path.join('pwmat', 'conf.lmp'), type_map = ['O', 'H'])
        tmp_system.to_pwmat_atomconfig('tmp.atom.config')
        self.system = dpdata.System()
        self.system.from_pwmat_atomconfig('tmp.atom.config')

class TestatomconfigDump1(unittest.TestCase, Testconfigoh):
    
    def setUp(self): 
        tmp_system = dpdata.System()
        tmp_system.from_pwmat_atomconfig(os.path.join('pwmat', 'atom.config.oh'))
        # tmp_system.from_lammps_lmp(os.path.join('poscars', 'conf.lmp'), type_map = ['O', 'H'])
        tmp_system.to_pwmat_atomconfig('tmp.atom.config')
        self.system = dpdata.System()
        self.system.from_pwmat_atomconfig('tmp.atom.config')

class TestatomconfigSkipZeroAtomNumb(unittest.TestCase) :
    def tearDown(self):
        if os.path.isfile('atom.config.tmp.1'):
            os.remove('atom.config.tmp.1')
        if os.path.isfile('atom.config.tmp.2'):
            os.remove('atom.config.tmp.2')

    def test_dump_pwmat_type_map(self):
        system0 = dpdata.System(os.path.join('pwmat', 'atom.config.oh'), fmt = 'pwmat/atom.config', type_map = ['H', 'O'])
        system0.to_pwmat_atomconfig('atom.config.tmp.1')
        system1 = dpdata.System(os.path.join('pwmat', 'atom.config.oh'), fmt = 'pwmat/atom.config', type_map = ['C', 'H', 'A', 'O', 'B'])
        system1.to_pwmat_atomconfig('atom.config.tmp.2')
        myfilecmp(self, 'atom.config.tmp.1', 'atom.config.tmp.2')


if __name__ == '__main__':
    unittest.main()
    
