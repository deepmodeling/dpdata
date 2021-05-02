import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompLabeledSys, IsPBC

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


if __name__ == '__main__':
    unittest.main()
