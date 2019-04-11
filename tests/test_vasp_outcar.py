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
        self.vir_places = 4

if __name__ == '__main__':
    unittest.main()
