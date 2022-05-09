import unittest
import numpy as np

from comp_sys import CompLabeledSys
from context import dpdata


@dpdata.driver.Driver.register("zero")
class ZeroDriver(dpdata.driver.Driver):
    def label(self, data):
        nframes = data['coords'].shape[0]
        natoms = data['coords'].shape[1]
        data['energies'] = np.zeros((nframes,))
        data['forces'] = np.zeros((nframes, natoms, 3))
        data['virials'] = np.zeros((nframes, 3, 3))
        return data


class TestPredict(unittest.TestCase, CompLabeledSys):
    def setUp (self) :
        ori_sys = dpdata.LabeledSystem('poscars/deepmd.h2o.md', 
                                        fmt = 'deepmd/raw', 
                                        type_map = ['O', 'H'])
        self.system_1 = ori_sys.predict(driver="zero")
        self.system_2 = dpdata.LabeledSystem('poscars/deepmd.h2o.md', 
                                             fmt = 'deepmd/raw', 
                                             type_map = ['O', 'H'])
        for pp in ('energies', 'forces', 'virials'):
            self.system_2.data[pp][:] = 0.

        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6
