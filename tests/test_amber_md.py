import os
import unittest
import shutil
from context import dpdata
from comp_sys import CompLabeledSys, IsPBC


class TestAmberMD(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp (self) :
        self.system_1 = dpdata.LabeledSystem('amber/02_Heat', fmt = 'amber/md')
        self.system_1.to('deepmd/npy','tmp.deepmd.npy')
        self.system_2 = dpdata.LabeledSystem('tmp.deepmd.npy', fmt = 'deepmd/npy')
        self.places = 5
        self.e_places = 4
        self.f_places = 6
        self.v_places = 6

    def tearDown(self) :
        if os.path.exists('tmp.deepmd.npy'):
            shutil.rmtree('tmp.deepmd.npy')

class TestAmberMDTarget(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        ll="amber/corr/low_level"
        ncfile="amber/corr/rc.nc"
        parmfile="amber/corr/qmmm.parm7"
        target = ":1"
        self.system_1 = dpdata.LabeledSystem(
            ll, nc_file=ncfile, parm7_file=parmfile, fmt='amber/md', use_element_symbols=target)
        self.system_2 = dpdata.LabeledSystem("amber/corr/dp_ll", fmt="deepmd/npy")

        self.places = 5
        self.e_places = 4
        self.f_places = 6
        self.v_places = 6

if __name__ == '__main__':
    unittest.main()