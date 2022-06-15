import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompLabeledSys, CompSys, IsNoPBC, IsPBC, MultiSystems

class TestDeepmdLoadDumpHDF5(unittest.TestCase, CompLabeledSys, IsPBC):
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


class TestDeepmdHDF5NoLabels(unittest.TestCase, CompSys, IsPBC) :
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


class TestHDF5Multi(unittest.TestCase, CompLabeledSys, MultiSystems, IsNoPBC):
    def setUp (self):
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

        system_1 = dpdata.LabeledSystem('gaussian/methane.gaussianlog', fmt='gaussian/log')
        system_2 = dpdata.LabeledSystem('gaussian/methane_reordered.gaussianlog', fmt='gaussian/log')
        system_3 = dpdata.LabeledSystem('gaussian/methane_sub.gaussianlog', fmt='gaussian/log')
        systems = dpdata.MultiSystems(system_1, system_2, system_3)
        systems.to_deepmd_hdf5("tmp.deepmd.hdf5")

        self.systems = dpdata.MultiSystems().from_deepmd_hdf5("tmp.deepmd.hdf5")
        self.system_names = ['C1H4', 'C1H3']
        self.system_sizes = {'C1H4':2, 'C1H3':1}
        self.atom_names = ['C', 'H']
        self.system_1 = self.systems['C1H3']
        self.system_2 = system_3

    def tearDown(self) :
        if os.path.exists('tmp.deepmd.hdf5'):
            os.remove('tmp.deepmd.hdf5')
