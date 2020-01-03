import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompSys, IsPBC
try:
   from pymatgen import  Structure
   exist_module=True
except:
   exist_module=False

@unittest.skipIf(not exist_module,"skip pymatgen")
class TestPymatgen(unittest.TestCase, CompSys, IsPBC):
    
    def setUp(self): 
        system_1 = dpdata.System()
        system_1.from_lammps_lmp(os.path.join('poscars', 'conf.lmp'), type_map = ['O', 'H'])
        system_1.to_pymatgen_structure()[0].to('poscar','tmp.POSCAR')
        self.system_1=system_1
        self.system_2=dpdata.System('tmp.POSCAR')
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

if __name__ == '__main__':
    unittest.main()
    
