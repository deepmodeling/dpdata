import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompSys
from comp_sys import CompLabeledSys
from comp_sys import IsPBC, IsNoPBC

class TestVaspXmlAppend(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp (self) :
        self.places = 6
        # rotated vasp computation, subject to numerical error
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6
        begin = 2
        end = 10
        step = 3
        self.system_1 = dpdata.LabeledSystem('poscars/vasprun.h2o.md.10.xml')
        self.system_2 = dpdata.LabeledSystem('poscars/vasprun.h2o.md.10.xml')
        self.system_1.append(self.system_2)
                
        self.system_1 = self.system_1.sub_system([0, 12, 4, 16, 8])
        self.system_2 = dpdata.LabeledSystem('poscars/vasprun.h2o.md.10.xml').sub_system(np.arange(0,10,2))


class TestDifferentOrderAppend(unittest.TestCase, CompLabeledSys, IsNoPBC):
    def setUp (self) :
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

        self.system_1 = dpdata.LabeledSystem('gaussian/methane.gaussianlog', fmt='gaussian/log')
        system_2 = dpdata.LabeledSystem('gaussian/methane_reordered.gaussianlog', fmt='gaussian/log')
        self.system_1.append(system_2)
                
        self.system_2 = self.system_1.sub_system([0, 0])

if __name__ == '__main__':
    unittest.main()
