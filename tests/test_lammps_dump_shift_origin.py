import unittest

from comp_sys import CompSys, IsPBC
from context import dpdata


class TestLammpsDumpShiftOrigin(unittest.TestCase, CompSys, IsPBC):
    def setUp(self):
        self.system_1 = dpdata.System("poscars/shift_origin.dump", fmt="lammps/dump")[0]
        self.system_2 = dpdata.System("poscars/shift_origin.dump", fmt="lammps/dump")[1]
        self.places = 6


if __name__ == "__main__":
    unittest.main()
