import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompSys
from comp_sys import CompLabeledSys
from comp_sys import IsPBC

class TestPWSCFTrajSkip(unittest.TestCase, CompSys, IsPBC):
    def setUp(self): 
        self.system_1 = dpdata.System(os.path.join('qe.traj', 'traj6'), 
                                      fmt = 'qe/cp/traj',
                                      begin = 1,
                                      step = 2)
        self.system_2 = dpdata.System(os.path.join('qe.traj', 'traj6'), 
                                      fmt = 'qe/cp/traj',
                                      begin = 0,
                                      step = 1) \
                              .sub_system(np.arange(1,6,2))
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4

class TestPWSCFLabeledTrajSkip(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self): 
        self.system_1 = dpdata.LabeledSystem(os.path.join('qe.traj', 'traj6'), 
                                             fmt = 'qe/cp/traj',
                                             begin = 1,
                                             step = 2)
        self.system_2 = dpdata.LabeledSystem(os.path.join('qe.traj', 'traj6'), 
                                             fmt = 'qe/cp/traj',
                                             begin = 0,
                                             step = 1) \
                              .sub_system(np.arange(1,6,2))
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4

    def test_cell(self):
        ref_cell = [5.359985500701728967e+00, 0, 0,
                    3.585941820098031974e-01, 5.317218997480877896e+00, 0,
                    7.606780476053129902e-01, 7.811107228901693622e-01, 5.715864930517207121e+00 ]
        ref_cell = 0.52917721067 * np.array(ref_cell).reshape(3,3)
        
        for ii in range(3) :
            for jj in range(3) :
                self.assertEqual(self.system_1.data['cells'][0][ii][jj], ref_cell[ii][jj])

        ref_cell = [5.308510801020571712e+00, 0, 0,
                    3.076052782312116429e-01, 5.279388982187173340e+00, 0,
                    4.321921336152507731e-01, 8.121110815096156399e-01, 5.301664983741235737e+00]
        ref_cell = 0.52917721067 * np.array(ref_cell).reshape(3,3)
        
        for ii in range(3) :
            for jj in range(3) :
                self.assertEqual(self.system_1.data['cells'][-1][ii][jj], ref_cell[ii][jj])
