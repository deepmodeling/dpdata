import unittest

from comp_sys import CompLabeledSys, IsPBC
from context import dpdata


class TestToList(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        system = dpdata.LabeledSystem("poscars/OUTCAR.h2o.md", fmt="vasp/outcar")
        self.system_1 = system.sub_system([2])
        self.system_2 = system.to_list()[2]
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4


if __name__ == "__main__":
    unittest.main()
