import unittest
from context import dpdata
from comp_sys import CompLabeledSys, IsPBC

class TestDeepmdLoadRaw(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp (self) :
        original_system = dpdata.LabeledSystem('poscars/OUTCAR.h2o.md', 
                                             fmt = 'vasp/outcar')
        original_system += original_system
        original_system += original_system
        original_system += original_system
        self.system_1 = dpdata.LabeledSystem()
        self.system_2 = original_system.copy()
        idx = self.system_2.shuffle()
        for ii in idx:
            self.system_1.append(original_system.sub_system(ii))

        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6
