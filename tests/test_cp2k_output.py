import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompLabeledSys

class TestCp2kNormalOutput(unittest.TestCase, CompLabeledSys):
    def setUp(self):
        self.system_1 = dpdata.LabeledSystem('cp2k/cp2k_normal_output/cp2k_output',fmt='cp2k/output')
        self.system_2 = dpdata.LabeledSystem('cp2k/cp2k_normal_output/deepmd', fmt='deepmd/npy')
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4
        
class TestCP2KDuplicateHeader(unittest.TestCase, CompLabeledSys):
    def setUp(self):
        self.system_1 = dpdata.LabeledSystem('cp2k/cp2k_duplicate_header/cp2k_output_duplicate_header',fmt='cp2k/output')
        self.system_2 = dpdata.LabeledSystem('cp2k/cp2k_duplicate_header/deepmd', fmt='deepmd/npy')                                      
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4

class TestCp2kReplaceElementOutput(unittest.TestCase, CompLabeledSys):
    def setUp(self):
        self.system_1 = dpdata.LabeledSystem('cp2k/cp2k_element_replace/cp2k_output_element_replace',fmt='cp2k/output')
        self.system_2 = dpdata.LabeledSystem('cp2k/cp2k_element_replace/deepmd', fmt='deepmd/npy')
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4

class TestNonCoveragedCP2KOutput(unittest.TestCase):
    def setUp (self) :
        self.system = dpdata.LabeledSystem('cp2k/cp2k_nocon_output',
                                           fmt = 'cp2k/output')
    def test_atom_types(self) :
        self.assertEqual(self.system.data['atom_types'], [])

    def test_cells(self) :
        self.assertEqual(self.system.data['cells'], [])

    def test_coords(self) :
        self.assertEqual(self.system.data['coords'], [])

    def test_energies(self) :
        self.assertEqual(self.system.data['energies'], [])

    def test_forces(self) :
        self.assertEqual(self.system.data['forces'], [])

    def test_virials(self) :
        self.assertFalse('virials' in self.system.data)


if __name__ == '__main__':
    unittest.main()

