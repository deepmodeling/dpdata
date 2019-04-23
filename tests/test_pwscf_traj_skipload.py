import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompSys
from comp_sys import CompLabeledSys

class TestPWSCFTrajSkip(unittest.TestCase, CompSys):
    def setUp(self): 
        self.system_1 = dpdata.System(os.path.join('pwscf.traj', 'traj6'), 
                                      fmt = 'pwscf/traj',
                                      begin = 1,
                                      step = 2)
        self.system_2 = dpdata.System(os.path.join('pwscf.traj', 'traj6'), 
                                      fmt = 'pwscf/traj',
                                      begin = 0,
                                      step = 1) \
                              .sub_system(np.arange(1,6,2))
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4

class TestPWSCFLabeledTrajSkip(unittest.TestCase, CompLabeledSys):
    def setUp(self): 
        self.system_1 = dpdata.LabeledSystem(os.path.join('pwscf.traj', 'traj6'), 
                                             fmt = 'pwscf/traj',
                                             begin = 1,
                                             step = 2)
        self.system_2 = dpdata.LabeledSystem(os.path.join('pwscf.traj', 'traj6'), 
                                             fmt = 'pwscf/traj',
                                             begin = 0,
                                             step = 1) \
                              .sub_system(np.arange(1,6,2))
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4

