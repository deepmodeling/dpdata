import os
import unittest

from context import dpdata


class TestGaussianGJF(unittest.TestCase):
    def setUp(self):
        self.system = dpdata.LabeledSystem("poscars/OUTCAR.h2o.md", fmt="vasp/outcar")

    def test_dump_gaussian_gjf(self):
        self.system.to_gaussian_gjf("tmp.gjf", keywords="force b3lyp/6-31g*")
        os.remove("tmp.gjf")
