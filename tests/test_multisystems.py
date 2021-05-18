import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompSys
from comp_sys import CompLabeledSys
from comp_sys import MultiSystems
from comp_sys import IsNoPBC
from itertools import permutations


class TestMultiSystems(unittest.TestCase, CompLabeledSys, MultiSystems, IsNoPBC):
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
    
        self.system_names = ['C1H4', 'C1H3']
        self.system_sizes = {'C1H4':2, 'C1H3':1}
        self.atom_names = ['C', 'H']


class TestMultiSystemsAdd(unittest.TestCase, CompLabeledSys, MultiSystems, IsNoPBC):
    def setUp(self):
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

        system_1 = dpdata.LabeledSystem('gaussian/methane.gaussianlog', fmt='gaussian/log')
        system_2 = dpdata.LabeledSystem('gaussian/methane_reordered.gaussianlog', fmt='gaussian/log')
        system_3 = dpdata.LabeledSystem('gaussian/methane_sub.gaussianlog', fmt='gaussian/log')
        system_4 = dpdata.LabeledSystem('gaussian/noncoveraged.gaussianlog', fmt='gaussian/log')

        self.systems = dpdata.MultiSystems(system_1)
        self.systems += system_2
        self.systems += system_3
        self.systems += system_4
        for s in self.systems:
            if s.formula == 'C1H3':
                self.system_1 = s
        self.system_2 = system_3
    
        self.system_names = ['C1H4', 'C1H3']
        self.system_sizes = {'C1H4':2, 'C1H3':1}
        self.atom_names = ['C', 'H']


class TestMultiSystemsSorted(unittest.TestCase, MultiSystems):
    def setUp(self):
        # CH4 and O2
        system_1 = dpdata.LabeledSystem('gaussian/methane.gaussianlog', fmt='gaussian/log')
        system_2 = dpdata.LabeledSystem('gaussian/oxygen.gaussianlog', fmt='gaussian/log')
        self.systems = dpdata.MultiSystems(system_1, system_2)

        self.system_names = ['C1H4O0', 'C0H0O2']
        self.system_sizes = {'C1H4O0':1, 'C0H0O2':1}
        self.atom_names = ['C', 'H', 'O']
        
class TestMultiDeepmdDumpRaw(unittest.TestCase, CompLabeledSys, IsNoPBC):
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

class TestMultiDeepmdDumpComp(unittest.TestCase, CompLabeledSys, IsNoPBC):
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

class TestTypeMap(unittest.TestCase):
    def setUp (self) :
        self.system_1 = dpdata.LabeledSystem('gaussian/methane.gaussianlog', fmt='gaussian/log')
        self.system_2 = dpdata.LabeledSystem('gaussian/methane_reordered.gaussianlog', fmt='gaussian/log')
        self.system_3 = dpdata.LabeledSystem('gaussian/methane_sub.gaussianlog', fmt='gaussian/log')
        self.system_4 = dpdata.LabeledSystem('gaussian/noncoveraged.gaussianlog', fmt='gaussian/log')

    def test_type_map(self):
        for type_map in permutations(['C', 'H', 'O', 'N'], 4):
            systems = dpdata.MultiSystems(self.system_1, self.system_2, self.system_3, self.system_4, type_map=type_map)
            self.assertEqual(type_map, systems.atom_names)


if __name__ == '__main__':
    unittest.main()
