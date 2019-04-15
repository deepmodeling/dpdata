import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompLabeledSys

class TestDeepmdLoadRaw(unittest.TestCase, CompLabeledSys):
    def setUp (self) :
        self.system_1 = dpdata.LabeledSystem('poscars/OUTCAR.h2o.md', 
                                             fmt = 'vasp/outcar')
        self.system_2 = dpdata.LabeledSystem('poscars/deepmd.h2o.md', 
                                             fmt = 'deepmd/raw', 
                                             type_map = ['O', 'H'])
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6


class TestDeepmdDumpRaw(unittest.TestCase, CompLabeledSys):
    def setUp (self) :
        self.system_1 = dpdata.LabeledSystem('poscars/OUTCAR.h2o.md', 
                                             fmt = 'vasp/outcar')
        self.system_1.to_deepmd_raw('tmp.deepmd')
        self.system_2 = dpdata.LabeledSystem('tmp.deepmd', type_map = ['O', 'H'])        
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6


if __name__ == '__main__':
    unittest.main()
