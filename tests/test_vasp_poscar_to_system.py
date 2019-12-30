import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompSys, IsPBC
from poscars.poscar_ref_oh import TestPOSCARoh

class TestPOSCARCart(unittest.TestCase, TestPOSCARoh):
    
    def setUp(self): 
        self.system = dpdata.System()
        self.system.from_vasp_poscar(os.path.join('poscars', 'POSCAR.oh.c'))

class TestPOSCARDirect(unittest.TestCase, TestPOSCARoh):
    
    def setUp(self): 
        self.system = dpdata.System()
        self.system.from_vasp_poscar(os.path.join('poscars', 'POSCAR.oh.d'))

class TestVaspPOSCARTypeMap(unittest.TestCase, CompSys, IsPBC):
    def setUp(self):
        sys0 = dpdata.System('poscars/POSCAR.oh.d', fmt =  'vasp/poscar')
        sys0.data['atom_names'] = ['A', 'H', 'B', 'O', 'D']
        sys0.data['atom_numbs'] = [  0,   1,   0,   1,   0]
        sys0.data['atom_types'] = np.array([  3, 1], dtype = int)
        sys1 = dpdata.System('poscars/POSCAR.oh.d', fmt =  'vasp/poscar', type_map = ['A', 'H', 'B', 'O', 'D'])
        self.system_1 = sys0
        self.system_2 = sys1
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6


if __name__ == '__main__':
    unittest.main()
