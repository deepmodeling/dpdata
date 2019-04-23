import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompLabeledSys

class TestVaspOUTCAR(unittest.TestCase, CompLabeledSys):
    def setUp (self) :
        self.system_1 = dpdata.LabeledSystem()
        self.system_1.from_vasp_xml('poscars/vasprun.h2o.md.xml')
        self.system_2 = dpdata.LabeledSystem()
        self.system_2.from_vasp_outcar('poscars/OUTCAR.h2o.md')
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4

class TestVaspOUTCARSkip(unittest.TestCase, CompLabeledSys):
    def setUp (self) :
        begin = 1
        step = 3
        end = 10
        self.system_1 = dpdata.LabeledSystem('poscars/OUTCAR.h2o.md.10', fmt = 'vasp/outcar', begin = begin, step = step)
        self.system_2 = dpdata.LabeledSystem('poscars/OUTCAR.h2o.md.10', fmt = 'vasp/outcar').sub_system(np.arange(begin, end, step))
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6


if __name__ == '__main__':
    unittest.main()
