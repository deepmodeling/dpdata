import os
import unittest
import shutil
from context import dpdata
from comp_sys import CompSys, CompLabeledSys, IsNoPBC

try:
    from dpdata import BondOrderSystem
except ImportError:
    skip_bond_order_system = True
else:
    skip_bond_order_system = False

class TestAmberSqmOut(unittest.TestCase, CompSys, IsNoPBC):
    def setUp (self) :
        self.system_1 = dpdata.System('amber/sqm_no_forces.out', fmt = 'sqm/out')
        self.system_1.to('deepmd/npy','tmp.sqm.noforces')
        self.system_2 = dpdata.System('tmp.sqm.noforces', fmt = 'deepmd/npy')
        self.places = 5
        self.e_places = 4
        self.f_places = 6
        self.v_places = 6

    def tearDown(self) :
        if os.path.exists('tmp.sqm.noforces'):
            shutil.rmtree('tmp.sqm.noforces')

class TestAmberSqmOutLabeled(unittest.TestCase, CompLabeledSys, IsNoPBC):
    def setUp(self) :
        self.system_1 = dpdata.LabeledSystem('amber/sqm_forces.out', fmt = 'sqm/out')
        self.system_1.to('deepmd/npy','tmp.sqm.forces')
        self.system_2 = dpdata.LabeledSystem('tmp.sqm.forces', fmt = 'deepmd/npy')
        self.places = 5
        self.e_places = 4
        self.f_places = 6
        self.v_places = 6

    def tearDown(self) :
        if os.path.exists('tmp.sqm.forces'):
            shutil.rmtree('tmp.sqm.forces')


class TestAmberSqmOutOpt(unittest.TestCase, CompLabeledSys, IsNoPBC):
    def setUp(self) :
        self.system_1 = dpdata.LabeledSystem('amber/sqm_opt.out', fmt = 'sqm/out')
        self.system_1.to('deepmd/npy','tmp.sqm.opt')
        self.system_2 = dpdata.LabeledSystem('tmp.sqm.opt', fmt = 'deepmd/npy')
        self.places = 5
        self.e_places = 4
        self.f_places = 6
        self.v_places = 6

    def tearDown(self) :
        if os.path.exists('tmp.sqm.opt'):
            shutil.rmtree('tmp.sqm.opt')


@unittest.skipIf(skip_bond_order_system, "dpdata does not have BondOrderSystem. One may install rdkit to fix.")
class TestAmberSqmIn(unittest.TestCase):
    def setUp(self):
        self.system = dpdata.BondOrderSystem("amber/methane.mol", fmt='mol', type_map=['H','C'])
        with open('amber/sqm.in', 'r') as f:
            self.sqm_in = f.read()
    
    def test_sqm_in(self):
        self.system.to("sqm/in", 'amber/sqm_test.in')
        with open('amber/sqm_test.in', 'r') as f:
            self.sqm_in_test = f.read()
        self.assertEqual(self.sqm_in, self.sqm_in_test)
    
    def tearDown(self):
        if os.path.isfile("amber/sqm_test.in"):
            os.remove("amber/sqm_test.in")

