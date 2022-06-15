import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompLabeledSys, IsPBC
from dpdata.utils import uniq_atom_names

class TestVaspOUTCAR(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp (self) :
        self.system_1 = dpdata.LabeledSystem()
        self.system_1.from_vasp_xml('poscars/vasprun.h2o.md.xml')
        self.system_2 = dpdata.LabeledSystem()
        self.system_2.from_vasp_outcar('poscars/OUTCAR.h2o.md')
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4

class TestVaspOUTCARTypeMap(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        sys0 = dpdata.LabeledSystem('poscars/OUTCAR.ch4.unconverged', fmt =  'vasp/outcar')
        sys0.data['atom_names'] = ['A', 'C', 'B', 'H', 'D']
        sys0.data['atom_numbs'] = [  0,   1,   0,   4,   0]
        sys0.data['atom_types'] = np.array([  3,   3,   3,   3,   1], dtype = int)
        sys1 = dpdata.LabeledSystem('poscars/OUTCAR.ch4.unconverged', fmt =  'vasp/outcar', type_map = ['A', 'C', 'B', 'H', 'D'])
        self.system_1 = sys0
        self.system_2 = sys1
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

class TestVaspOUTCARSkip(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp (self) :
        begin = 1
        step = 3
        end = 10
        self.system_1 = dpdata.LabeledSystem('poscars/OUTCAR.h2o.md.10', fmt = 'vasp/outcar', begin = begin, step = step)
        self.system_2 = dpdata.LabeledSystem('poscars/OUTCAR.h2o.md.10', fmt = 'vasp/outcar').sub_system(np.arange(begin, end, step))
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6


class TestVaspOUTCARVdw(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp (self) :
        self.system_1 = dpdata.LabeledSystem('poscars/OUTCAR.Ge.vdw', fmt = 'vasp/outcar')
        self.system_2 = dpdata.LabeledSystem()
        self.system_2.from_vasp_xml('poscars/vasprun.Ge.vdw.xml')
        self.places = 5
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6


class TestDuplicatedAtomNames(unittest.TestCase):
    def test(self):
        system = dpdata.LabeledSystem('poscars/6362_OUTCAR', fmt = 'vasp/outcar')
        expected_types = [0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1]
        self.assertEqual(list(system['atom_types']), expected_types)
        self.assertEqual(system['atom_names'], ['B', 'O'])
        self.assertEqual(system['atom_numbs'], [8, 6])

    def test_type_map(self):
        system = dpdata.LabeledSystem('poscars/6362_OUTCAR', fmt = 'vasp/outcar', type_map = ['O', 'B'])
        expected_types = [1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0]
        self.assertEqual(list(system['atom_types']), expected_types)
        self.assertEqual(system['atom_names'], ['O', 'B'])
        self.assertEqual(system['atom_numbs'], [6, 8])


class TestUniqAtomNames(unittest.TestCase):
    def test(self):
        data = {}
        data['atom_names'] = ['O', 'H', 'O', 'H']
        data['atom_types'] = np.array([0, 1, 2, 3, 3, 2, 1], dtype=int)
        
        data = uniq_atom_names(data)
        self.assertEqual(list(data['atom_types']),
                         [0, 1, 0, 1, 1, 0, 1])
        self.assertEqual(list(data['atom_names']),
                         ['O', 'H'])
        self.assertEqual(list(data['atom_numbs']),
                         [3, 4])

class TestVaspOUTCARML(unittest.TestCase):
    def test(self):
        system1 = dpdata.LabeledSystem('poscars/OUTCAR.ch4.ml', fmt = 'vasp/outcar',ml=True)
        system2 = dpdata.LabeledSystem('poscars/OUTCAR.ch4.ml', fmt = 'vasp/outcar',ml=False)
        expected_types = [0, 0, 0, 0, 1]
        self.assertEqual(list(system1['atom_types']), expected_types)
        self.assertEqual(system1['atom_names'], ['H', 'C'])
        self.assertEqual(len(system1['energies']), 10)
        self.assertEqual(list(system2['atom_types']), expected_types)
        self.assertEqual(system2['atom_names'], ['H', 'C'])
        self.assertEqual(len(system2['energies']), 4)


if __name__ == '__main__':
    unittest.main()
