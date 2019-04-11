import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompSys

class TestVaspXml(unittest.TestCase, CompSys):
    def setUp (self) :
        self.places = 6
        xml_sys = dpdata.LabeledSystem()
        xml_sys.from_vasp_xml('poscars/vasprun.h2o.md.xml')
        # init_sys = dpdata.System()
        # init_sys.from_vasp_poscar('poscars/POSCAR.h2o.md')
        finl_sys = dpdata.System()
        finl_sys.from_vasp_poscar('poscars/CONTCAR.h2o.md')
        self.system_1 = finl_sys
        self.system_2 = xml_sys.sub_system([-1])

if __name__ == '__main__':
    unittest.main()
