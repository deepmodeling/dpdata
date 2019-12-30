import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompSys, IsPBC

class TestLmpDumpSkip(unittest.TestCase, CompSys, IsPBC):
    
    def setUp(self): 
        self.system_1 = dpdata.System(os.path.join('poscars', 'conf.5.dump'), 
                                      type_map = ['O', 'H'], 
                                      begin = 1,
                                      step = 2)
        self.system_2 = dpdata.System(os.path.join('poscars', 'conf.5.dump'), 
                                      type_map = ['O', 'H'], 
                                      begin = 0,
                                      step = 1) \
                              .sub_system(np.arange(1,5,2))
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4
