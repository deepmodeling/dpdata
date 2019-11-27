import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompLabeledSys

class TestGapxyz(unittest.TestCase, CompLabeledSys):
    def setUp (self) :
        self.multi_systems = dpdata.MultiSystems(file_name='xyz/xyz_unittest.xyz', fmt='gap/xyz')
        self.system_1 = self.multi_systems.systems['B1C9']
        self.system_2 = dpdata.LabeledSystem('xyz/B1C9', fmt='deepmd')
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4

class TestGapxyz2(unittest.TestCase, CompLabeledSys):
    def setUp (self) :
        self.system_temp0 = dpdata.MultiSystems(file_name='xyz/xyz_unittest.xyz', fmt='gap/xyz')
        self.system_1 = self.system_temp0.systems['B5C7']
        self.system_temp1 = dpdata.LabeledSystem('xyz/B1C9', fmt='deepmd')
        self.system_temp2 = dpdata.LabeledSystem('xyz/B5C7', fmt='deepmd')
        self.system_temp3 = dpdata.MultiSystems(self.system_temp1, self.system_temp2)
        self.system_2 = self.system_temp3.systems['B5C7']
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4

if __name__ == '__main__':
    unittest.main()
