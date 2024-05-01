import subprocess as sp
import sys
import unittest

from context import dpdata
from poscars.poscar_ref_oh import TestPOSCARoh


class TestCli(unittest.TestCase, TestPOSCARoh):
    @classmethod
    def setUpClass(cls) -> None:
        sp.check_output(
            [
                "dpdata",
                "poscars/conf.lmp",
                "--type-map",
                "O",
                "H",
                "-olammps/lmp",
                "-O",
                "tmp.lmp",
                "--no-labeled",
            ]
        )
        cls.system = dpdata.System("tmp.lmp", fmt="lammps/lmp", type_map=["O", "H"])

    @classmethod
    def tearDownClass(cls) -> None:
        cls.system = None


class TestClassScript(unittest.TestCase):
    def test_class_script(self):
        expected_version = dpdata.__version__
        output = sp.check_output([sys.executable, "-m", "dpdata", "--version"]).decode(
            "ascii"
        )
        assert output.splitlines()[0] == f"dpdata v{expected_version}"
