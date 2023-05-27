import os
import unittest

from comp_sys import CompSys
from context import dpdata


class TestGaussianGJF(unittest.TestCase):
    def setUp(self):
        self.system = dpdata.LabeledSystem("poscars/OUTCAR.h2o.md", fmt="vasp/outcar")

    def test_dump_gaussian_gjf(self):
        self.system.to_gaussian_gjf("tmp.gjf", keywords="force b3lyp/6-31g*")
        os.remove("tmp.gjf")


class TestGaussianGJFComp(unittest.TestCase, CompSys):
    def setUp(self):
        self.system_1 = dpdata.LabeledSystem(
            "poscars/OUTCAR.h2o.md", fmt="vasp/outcar"
        )[0]
        self.system_1.to_gaussian_gjf("tmp.gjf", keywords="force b3lyp/6-31g*")
        self.system_2 = dpdata.System(
            "tmp.gjf", fmt="gaussian/gjf", type_map=self.system_1.get_atom_names()
        )
        os.remove("tmp.gjf")
        self.places = 6
