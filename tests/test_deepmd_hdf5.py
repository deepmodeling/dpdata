import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompLabeledSys, CompSys, IsPBC

class TestDeepmdLoadDumpComp(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp (self) :
        self.system_1 = dpdata.LabeledSystem('poscars/OUTCAR.h2o.md', 
                                             fmt = 'vasp/outcar')
        self.system_1.to_deepmd_hdf5('tmp.deepmd.hdf5', 
                                    prec = np.float64, 
                                    set_size = 2)        

        self.system_2 = dpdata.LabeledSystem('tmp.deepmd.hdf5', 
                                             fmt = 'deepmd/hdf5',
                                             type_map = ['O', 'H'])
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

    def tearDown(self) :
        if os.path.exists('tmp.deepmd.hdf5'):
            os.remove('tmp.deepmd.hdf5')


class TestDeepmdCompNoLabels(unittest.TestCase, CompSys, IsPBC) :
    def setUp (self) :
        self.system_1 = dpdata.System('poscars/POSCAR.h2o.md',
                                      fmt = 'vasp/poscar')
        self.system_1.to_deepmd_hdf5('tmp.deepmd.hdf5', 
                                    prec = np.float64, 
                                    set_size = 2)        
        self.system_2 = dpdata.System('tmp.deepmd.hdf5', 
                                      fmt = 'deepmd/hdf5',
                                      type_map = ['O', 'H'])
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

    def tearDown(self) :
        if os.path.exists('tmp.deepmd.hdf5'):
            os.remove('tmp.deepmd.hdf5')
