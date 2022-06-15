import unittest
import numpy as np

from comp_sys import CompLabeledSys, IsPBC
from context import dpdata
try:
    import ase
except ModuleNotFoundError:
    skip_ase = True
else:
    skip_ase = False


@dpdata.driver.Driver.register("zero")
class ZeroDriver(dpdata.driver.Driver):
    def label(self, data):
        nframes = data['coords'].shape[0]
        natoms = data['coords'].shape[1]
        data['energies'] = np.zeros((nframes,))
        data['forces'] = np.zeros((nframes, natoms, 3))
        data['virials'] = np.zeros((nframes, 3, 3))
        return data


@dpdata.driver.Driver.register("one")
class OneDriver(dpdata.driver.Driver):
    def label(self, data):
        nframes = data['coords'].shape[0]
        natoms = data['coords'].shape[1]
        data['energies'] = np.ones((nframes,))
        data['forces'] = np.ones((nframes, natoms, 3))
        data['virials'] = np.ones((nframes, 3, 3))
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


class TestHybridDriver(unittest.TestCase, CompLabeledSys):
    """Test HybridDriver."""
    def setUp(self) :
        ori_sys = dpdata.LabeledSystem('poscars/deepmd.h2o.md', 
                                        fmt = 'deepmd/raw', 
                                        type_map = ['O', 'H'])
        self.system_1 = ori_sys.predict([
                {"type": "one"},
                {"type": "one"},
                {"type": "one"},
                {"type": "zero"},
            ],
            driver="hybrid")
        # sum is 3
        self.system_2 = dpdata.LabeledSystem('poscars/deepmd.h2o.md', 
                                             fmt = 'deepmd/raw', 
                                             type_map = ['O', 'H'])
        for pp in ('energies', 'forces'):
            self.system_2.data[pp][:] = 3.

        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6


@unittest.skipIf(skip_ase,"skip ase related test. install ase to fix")
class TestASEtraj1(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp (self) :
        ori_sys = dpdata.LabeledSystem('poscars/deepmd.h2o.md', 
                                        fmt = 'deepmd/raw', 
                                        type_map = ['O', 'H'])
        one_driver = OneDriver()
        self.system_1 = ori_sys.predict(driver=one_driver)
        self.system_2 = ori_sys.predict(one_driver.ase_calculator, driver="ase")
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4
