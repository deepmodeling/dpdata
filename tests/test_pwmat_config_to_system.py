import os
import numpy as np
import unittest
import dpdata
from pwmat.config_ref_ch4 import Testconfigch4

class Testconfig(unittest.TestCase, Testconfigch4):
    
    def setUp(self): 
        self.system = dpdata.System()
        self.system.from_pwmat_atomconfig(os.path.join('pwmat', 'atom.config'))
class TestpwmatconfigTypeMap(unittest.TestCase):
    def setUp(self):
        sys0 = dpdata.System('pwmat/atom.config', fmt =  'atom.config')
        sys0.data['atom_names'] = ['A', 'H', 'B', 'C', 'D']
        sys0.data['atom_numbs'] = [  0,   1,   0,   1,   0]
        sys0.data['atom_types'] = np.array([  0, 0, 0, 1], dtype = int)
        sys1 = dpdata.System('pwmat/atom.config', fmt =  'pwmat/atom.config', type_map = ['A', 'H', 'B', 'C', 'D'])
        self.system_1 = sys0
        self.system_2 = sys1
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6


if __name__ == '__main__':
    unittest.main()
