import unittest
from context import dpdata
from poscars.poscar_ref_oh import TestPOSCARoh
import subprocess as sp


class TestCli(unittest.TestCase, TestPOSCARoh):
    
    @classmethod
    def setUpClass(cls) -> None:
        sp.check_output(["dpdata", "poscars/conf.lmp", "--type-map", "O", "H", "-olammps/lmp", "-O", "tmp.lmp", "--no-labeled"])
        cls.system = dpdata.System('tmp.lmp', fmt='lammps/lmp',
                             type_map = ['O', 'H'])

    @classmethod
    def tearDownClass(cls) -> None:
        cls.system = None
