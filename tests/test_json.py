import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompLabeledSys, IsPBC

class TestJsonLoad(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp (self) :
        self.system_1 = dpdata.LabeledSystem('poscars/OUTCAR.h2o.md',
                                             fmt = 'vasp/outcar')
        self.system_2 = dpdata.LabeledSystem.load('poscars/h2o.md.json')
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4

class TestAsDict(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp (self) :
        self.system_1 = dpdata.LabeledSystem('poscars/OUTCAR.h2o.md',
                                             fmt = 'vasp/outcar')
        self.system_2 = dpdata.LabeledSystem.from_dict(self.system_1.as_dict())
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4

if __name__ == '__main__':
    unittest.main()
