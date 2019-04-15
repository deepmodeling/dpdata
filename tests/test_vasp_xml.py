import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompSys
from comp_sys import CompLabeledSys

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


class TestVaspXmlRotSys(unittest.TestCase, CompLabeledSys):
    def setUp (self) :
        self.places = 4
        # rotated vasp computation, subject to numerical error
        self.e_places = 3
        self.f_places = 2
        self.v_places = 1
        self.system_1 = dpdata.LabeledSystem('poscars/vasprun.h2o.md.tribox.xml')
        self.system_2 = dpdata.LabeledSystem('poscars/vasprun.h2o.md.tribox.lower.xml')


if __name__ == '__main__':
    unittest.main()
