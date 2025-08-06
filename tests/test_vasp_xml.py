from __future__ import annotations

import unittest

import numpy as np
from comp_sys import CompLabeledSys, CompSys, IsPBC
from context import dpdata


class TestVaspXml(unittest.TestCase, CompSys, IsPBC):
    def setUp(self):
        self.places = 6
        xml_sys = dpdata.LabeledSystem()
        xml_sys.from_vasp_xml("poscars/vasprun.h2o.md.xml")
        # init_sys = dpdata.System()
        # init_sys.from_vasp_poscar('poscars/POSCAR.h2o.md')
        finl_sys = dpdata.System()
        finl_sys.from_vasp_poscar("poscars/CONTCAR.h2o.md")
        self.system_1 = finl_sys
        self.system_2 = xml_sys.sub_system([-1])


class TestVaspXmlConvTrue(unittest.TestCase, CompSys, IsPBC):
    def setUp(self):
        self.system_1 = dpdata.LabeledSystem()
        # the first frame is not converged.
        self.system_1.from_vasp_xml(
            "poscars/vasprun.h2o.md.conv.xml", convergence_check=True
        )
        self.system_2 = dpdata.LabeledSystem()
        self.system_2.from_vasp_xml("poscars/vasprun.h2o.md.xml")
        # check if frames 1:3 match
        self.system_2 = self.system_2[1:]
        self.places = 6


class TestVaspXmlConvFalse(unittest.TestCase, CompSys, IsPBC):
    def setUp(self):
        self.system_1 = dpdata.LabeledSystem()
        self.system_1.from_vasp_xml(
            "poscars/vasprun.h2o.md.conv.xml", convergence_check=False
        )
        self.system_2 = dpdata.LabeledSystem()
        self.system_2.from_vasp_xml("poscars/vasprun.h2o.md.xml")
        self.places = 6


class TestVaspXmlDup(unittest.TestCase, CompSys, IsPBC):
    def setUp(self):
        self.places = 6
        xml_sys = dpdata.LabeledSystem()
        xml_sys.from_vasp_xml("poscars/vasprun.h2o.md.duplicate.xml")
        finl_sys = dpdata.System()
        finl_sys.from_vasp_poscar("poscars/CONTCAR.h2o.md")
        self.system_1 = finl_sys
        self.system_2 = xml_sys.sub_system([-1])


class TestVaspXmlRotSys(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        self.places = 4
        # rotated vasp computation, subject to numerical error
        self.e_places = 3
        self.f_places = 2
        self.v_places = 1
        self.system_1 = dpdata.LabeledSystem("poscars/vasprun.h2o.md.tribox.xml")
        self.system_2 = dpdata.LabeledSystem("poscars/vasprun.h2o.md.tribox.lower.xml")


class TestVaspXmlSkip(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        self.places = 6
        # rotated vasp computation, subject to numerical error
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6
        begin = 2
        end = 10
        step = 3
        self.system_1 = dpdata.LabeledSystem(
            "poscars/vasprun.h2o.md.10.xml", begin=begin, step=step
        )
        self.system_2 = dpdata.LabeledSystem(
            "poscars/vasprun.h2o.md.10.xml"
        ).sub_system(np.arange(2, 10, 3))


class TestVaspXmlNoVirial(unittest.TestCase, CompSys, IsPBC):
    def setUp(self):
        self.places = 6
        xml_sys = dpdata.LabeledSystem()
        xml_sys.from_vasp_xml("poscars/vasprun.h2o.md.novirial.xml")
        finl_sys = dpdata.System()
        finl_sys.from_vasp_poscar("poscars/CONTCAR.h2o.md")
        self.system_1 = finl_sys
        self.system_2 = xml_sys.sub_system([-1])


class TestVaspOUTCARNWRITE0(unittest.TestCase):
    def test(self):
        # all frames are written to vasprun.xml that have forces are read
        # even though the nwrite parameter is set to 0
        ss = dpdata.LabeledSystem("poscars/Ti-aimd-nwrite0/vasprun.xml")
        self.assertEqual(ss.get_nframes(), 10)


class TestVaspAtomNamesV6(unittest.TestCase):
    def test(self):
        ss = dpdata.LabeledSystem("poscars/Ti-O-Ti-v6/vasprun.xml")
        self.assertEqual(ss.get_atom_names(), ["Ti", "O"])
        self.assertEqual(ss.get_atom_numbs(), [6, 2])
        np.testing.assert_equal(ss.get_atom_types(), [0, 0, 0, 1, 1, 0, 0, 0])


if __name__ == "__main__":
    unittest.main()
