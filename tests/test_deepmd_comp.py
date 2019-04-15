import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompLabeledSys

class TestDeepmdLoadDumpComp(unittest.TestCase, CompLabeledSys):
    def setUp (self) :
        self.system_1 = dpdata.LabeledSystem('poscars/OUTCAR.h2o.md', 
                                             fmt = 'vasp/outcar')
        self.system_1.to_deepmd_npy('tmp.deepmd.npy', 
                                    prec = np.float64, 
                                    set_size = 2)        

        self.system_2 = dpdata.LabeledSystem('tmp.deepmd.npy', 
                                             fmt = 'deepmd/npy',
                                             type_map = ['O', 'H'])
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6


if __name__ == '__main__':
    unittest.main()
