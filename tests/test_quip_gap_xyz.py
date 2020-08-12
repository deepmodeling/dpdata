import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompLabeledSys, IsPBC

class TestQuipGapxyz1(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp (self) :
        self.multi_systems = dpdata.MultiSystems.from_file('xyz/xyz_unittest.xyz','quip/gap/xyz')
        self.system_1 = self.multi_systems.systems['B1C9']
        self.system_2 = dpdata.LabeledSystem('xyz/B1C9', fmt='deepmd')
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4

class TestQuipGapxyz2(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp (self) :
        self.system_temp0 = dpdata.MultiSystems.from_file(file_name='xyz/xyz_unittest.xyz', fmt='quip/gap/xyz')
        self.system_1 = self.system_temp0.systems['B5C7'] # .sort_atom_types()
        self.system_temp1 = dpdata.LabeledSystem('xyz/B1C9', fmt='deepmd')
        self.system_temp2 = dpdata.LabeledSystem('xyz/B5C7', fmt='deepmd')
        self.system_temp3 = dpdata.MultiSystems(self.system_temp2, self.system_temp1)
        self.system_2 = self.system_temp3.systems['B5C7']
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4

class TestQuipGapxyzsort1(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp (self) :
        self.multi_systems_1 = dpdata.MultiSystems.from_file('xyz/xyz_unittest.sort.xyz','quip/gap/xyz')
        self.system_1 = self.multi_systems_1.systems['B5C7']
        self.system_1.sort_atom_types()
        self.multi_systems_2 = dpdata.MultiSystems.from_file('xyz/xyz_unittest.xyz','quip/gap/xyz')
        self.system_2 = self.multi_systems_2.systems['B5C7']
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4

class TestQuipGapxyzsort2(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp (self) :
        self.multi_systems_1 = dpdata.MultiSystems.from_file('xyz/xyz_unittest.sort.xyz','quip/gap/xyz')
        self.system_1 = self.multi_systems_1.systems['B1C9']
        self.system_1.sort_atom_types()
        self.multi_systems_2 = dpdata.MultiSystems.from_file('xyz/xyz_unittest.xyz','quip/gap/xyz')
        self.system_2 = self.multi_systems_2.systems['B1C9']
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4

class TestQuipGapxyzfield(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp (self) :
        self.multi_systems_1 = dpdata.MultiSystems.from_file('xyz/xyz_unittest.field.xyz','quip/gap/xyz')
        self.system_1 = self.multi_systems_1.systems['B1C9']
        self.system_1.sort_atom_types()
        self.multi_systems_2 = dpdata.MultiSystems.from_file('xyz/xyz_unittest.xyz','quip/gap/xyz')
        self.system_2 = self.multi_systems_2.systems['B1C9']
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4

class TestQuipGapxyzfield2(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp (self) :
        self.multi_systems_1 = dpdata.MultiSystems.from_file('xyz/xyz_unittest.field.xyz','quip/gap/xyz')
        self.system_1 = self.multi_systems_1.systems['B5C7']
        self.system_1.sort_atom_types()
        self.multi_systems_2 = dpdata.MultiSystems.from_file('xyz/xyz_unittest.xyz','quip/gap/xyz')
        self.system_2 = self.multi_systems_2.systems['B5C7']
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4

class TestQuipGapxyzNoVirials(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        self.multi_systems_1 = dpdata.MultiSystems.from_file('xyz/xyz_B5C7_novirials.xyz', fmt='quip/gap/xyz')
        self.system_1 = self.multi_systems_1.systems['B5C7']
        self.system_1.sort_atom_types()
        self.system_2 = dpdata.LabeledSystem('xyz/B5C7_novirials', fmt='deepmd/raw')
        self.places = 6
        self.e_places = 6
        self.f_places = 6



if __name__ == '__main__':
    unittest.main()
