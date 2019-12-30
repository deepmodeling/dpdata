import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompSys
from comp_sys import CompLabeledSys
from comp_sys import IsPBC

class TestVaspXml(unittest.TestCase, CompSys, IsPBC):
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


class TestVaspXmlRotSys(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp (self) :
        self.places = 4
        # rotated vasp computation, subject to numerical error
        self.e_places = 3
        self.f_places = 2
        self.v_places = 1
        self.system_1 = dpdata.LabeledSystem('poscars/vasprun.h2o.md.tribox.xml')
        self.system_2 = dpdata.LabeledSystem('poscars/vasprun.h2o.md.tribox.lower.xml')


class TestVaspXmlSkip(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp (self) :
        self.places = 6
        # rotated vasp computation, subject to numerical error
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6
        begin = 2
        end = 10
        step = 3
        self.system_1 = dpdata.LabeledSystem('poscars/vasprun.h2o.md.10.xml', begin = begin, step = step)
        self.system_2 = dpdata.LabeledSystem('poscars/vasprun.h2o.md.10.xml').sub_system(np.arange(2,10,3))


if __name__ == '__main__':
    unittest.main()
