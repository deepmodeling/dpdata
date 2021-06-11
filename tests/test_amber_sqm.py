import os
import unittest
import shutil
from context import dpdata
from comp_sys import CompSys, IsNoPBC


class TestAmberSqmOut(unittest.TestCase, CompSys, IsNoPBC):
    def setUp (self) :
        self.system_1 = dpdata.System('amber/sqm.out', fmt = 'sqm/out')
        self.system_1.to('deepmd/npy','tmp.deepmd.npy')
        self.system_2 = dpdata.System('tmp.deepmd.npy', fmt = 'deepmd/npy')
        self.places = 5
        self.e_places = 4
        self.f_places = 6
        self.v_places = 6

    def tearDown(self) :
        if os.path.exists('tmp.deepmd.npy'):
            shutil.rmtree('tmp.deepmd.npy')

