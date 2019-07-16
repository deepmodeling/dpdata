import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompSys
from comp_sys import CompLabeledSys

class TestMultiSystems(unittest.TestCase, CompLabeledSys):
    def setUp(self):
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

        system_1 = dpdata.LabeledSystem('gaussian/methane.gaussianlog', fmt='gaussian/log')
        system_2 = dpdata.LabeledSystem('gaussian/methane_reordered.gaussianlog', fmt='gaussian/log')
        system_3 = dpdata.LabeledSystem('gaussian/methane_sub.gaussianlog', fmt='gaussian/log')
        system_4 = dpdata.LabeledSystem('gaussian/noncoveraged.gaussianlog', fmt='gaussian/log')

        self.systems = dpdata.MultiSystems(system_1, system_3, system_4)
        self.systems.append(system_2)
        self.system_1 = self.systems['C1H3']
        self.system_2 = system_3
    
    def test_systems_name(self):
        system_names = ['C1H4', 'C1H3']
        self.assertEqual(set(self.systems.systems), set(system_names))
    
    def test_systems_size(self):
        system_sizes = {'C1H4':2, 'C1H3':1}
        for name, size in system_sizes.items():
            self.assertEqual(self.systems[name].get_nframes(), size)

class TestMultiDeepmdDumpRaw(unittest.TestCase, CompLabeledSys):
    def setUp (self) :
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

        system_1 = dpdata.LabeledSystem('gaussian/methane.gaussianlog', fmt='gaussian/log')
        system_2 = dpdata.LabeledSystem('gaussian/methane_reordered.gaussianlog', fmt='gaussian/log')
        system_3 = dpdata.LabeledSystem('gaussian/methane_sub.gaussianlog', fmt='gaussian/log')
        system_4 = dpdata.LabeledSystem('gaussian/noncoveraged.gaussianlog', fmt='gaussian/log')

        systems = dpdata.MultiSystems(system_1, system_2, system_3, system_4)
        path = "tmp.deepmd.multi"
        systems.to_deepmd_raw(path)
        self.system_1 = dpdata.LabeledSystem(os.path.join(path, 'C1H3'), fmt='deepmd/raw', type_map = ['C', 'H'])
        self.system_2 = system_3

class TestMultiDeepmdDumpComp(unittest.TestCase, CompLabeledSys):
    def setUp (self) :
        self.places = 6
        self.e_places = 4
        self.f_places = 6
        self.v_places = 6

        system_1 = dpdata.LabeledSystem('gaussian/methane.gaussianlog', fmt='gaussian/log')
        system_2 = dpdata.LabeledSystem('gaussian/methane_reordered.gaussianlog', fmt='gaussian/log')
        system_3 = dpdata.LabeledSystem('gaussian/methane_sub.gaussianlog', fmt='gaussian/log')
        system_4 = dpdata.LabeledSystem('gaussian/noncoveraged.gaussianlog', fmt='gaussian/log')

        systems = dpdata.MultiSystems(system_1, system_2, system_3, system_4)
        path = "tmp.deepmd.npy.multi"
        systems.to_deepmd_npy(path)
        self.system_1 = dpdata.LabeledSystem(os.path.join(path, 'C1H3'), fmt='deepmd/npy', type_map = ['C', 'H'])
        self.system_2 = system_3

if __name__ == '__main__':
    unittest.main()
