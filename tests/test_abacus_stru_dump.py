from __future__ import annotations

import os
import shutil
import unittest

from context import dpdata
from test_vasp_poscar_dump import myfilecmp

from dpdata.abacus.scf import parse_stru_pos


class TestStruDump(unittest.TestCase):
    def setUp(self):
        self.system_ch4 = dpdata.System("abacus.scf/STRU.ch4", fmt="stru")

    def tearDown(self):
        if os.path.isfile("STRU_tmp"):
            os.remove("STRU_tmp")

    def test_dump_stru(self):
        self.system_ch4.to(
            "stru",
            "STRU_tmp",
            mass=[12, 1],
            pp_file=["C.upf", "H.upf"],
            numerical_orbital=["C.orb", "H.orb"],
            numerical_descriptor="jle.orb",
        )
        myfilecmp(self, "abacus.scf/stru_test", "STRU_tmp")

    def test_dumpStruLinkFile(self):
        os.makedirs("abacus.scf/tmp", exist_ok=True)
        self.system_ch4.to(
            "stru",
            "abacus.scf/tmp/STRU_tmp",
            mass=[12, 1],
            pp_file=["abacus.scf/C.upf", "abacus.scf/H.upf"],
            numerical_orbital=["abacus.scf/C.orb", "abacus.scf/H.orb"],
            numerical_descriptor="abacus.scf/jle.orb",
            link_file=True,
        )
        myfilecmp(self, "abacus.scf/stru_test", "abacus.scf/tmp/STRU_tmp")

        self.assertTrue(os.path.islink("abacus.scf/tmp/C.upf"))
        self.assertTrue(os.path.islink("abacus.scf/tmp/H.upf"))
        self.assertTrue(os.path.islink("abacus.scf/tmp/C.orb"))
        self.assertTrue(os.path.islink("abacus.scf/tmp/H.orb"))
        self.assertTrue(os.path.islink("abacus.scf/tmp/jle.orb"))

        if os.path.isdir("abacus.scf/tmp"):
            shutil.rmtree("abacus.scf/tmp")
    
    def test_dump_stru_pporb_mismatch(self):
        with self.assertRaises(KeyError, msg="pp_file is a dict and lack of pp for H"):
            self.system_ch4.to(
                "stru",
                "STRU_tmp",
                mass=[12, 1],
                pp_file={"C": "C.upf", "O": "O.upf"},
                numerical_orbital={"C": "C.orb", "H": "H.orb"},
            )
            
        with self.assertRaises(ValueError, msg="pp_file is a list and lack of pp for H"):
            self.system_ch4.to(
                "stru",
                "STRU_tmp",
                mass=[12, 1],
                pp_file=["C.upf"],
                numerical_orbital={"C": "C.orb", "H": "H.orb"},
            )
    
        with self.assertRaises(KeyError, msg="numerical_orbital is a dict and lack of orbital for H"):
            self.system_ch4.to(
                "stru",
                "STRU_tmp",
                mass=[12, 1],
                pp_file={"C": "C.upf", "H": "H.upf"},
                numerical_orbital={"C": "C.orb", "O": "O.orb"},
            )
            
        with self.assertRaises(ValueError, msg="numerical_orbital is a list and lack of orbital for H"):
            self.system_ch4.to(
                "stru",
                "STRU_tmp",
                mass=[12, 1],
                pp_file=["C.upf", "H.upf"],
                numerical_orbital=["C.orb"],
            )

    def test_dump_spinconstrain(self):
        self.system_ch4.to(
            "stru",
            "STRU_tmp",
            mass=[12, 1],
            pp_file={"C": "C.upf", "H": "H.upf"},
            numerical_orbital={"C": "C.orb", "H": "H.orb"},
            mag=[4, [1, 1, 1], 1, 1, 1],
            sc=[True, True, [True, False, True], False, True],
            move=[1, 1, 1, 0, 0],
            angle1=[None, None, 100, 90, 80],
            angle2=[100, 90, 80, 70, None],
            lambda_=[[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9], None, None],
        )

        assert os.path.isfile("STRU_tmp")
        with open("STRU_tmp") as f:
            c = f.read()

        with open("abacus.scf/stru.ref") as f:
            stru_ref = f.read()
        assert c == stru_ref

    def test_dump_spin(self):
        sys_tmp = dpdata.System("abacus.scf/stru.ref", fmt="stru")
        sys_tmp.data["spins"] = [
            [[1, 2, 3], [4, 5, 6], [1, 1, 1], [2, 2, 2], [3, 3, 3]]
        ]
        sys_tmp.to(
            "stru",
            "STRU_tmp",
            mass=[12, 1],
            pp_file=["C.upf", "H.upf"],
            numerical_orbital=["C.orb", "H.orb"],
        )
        assert os.path.isfile("STRU_tmp")
        with open("STRU_tmp") as f:
            c = f.read()
        stru_ref = """C
0.0
1
5.192682633809 4.557725978258 4.436846615358 1 1 1 mag 1.000000000000 2.000000000000 3.000000000000
H
0.0
4
5.416431453540 4.011298860305 3.511161492417 1 1 1 mag 4.000000000000 5.000000000000 6.000000000000
4.131588222365 4.706745191323 4.431136645083 1 1 1 mag 1.000000000000 1.000000000000 1.000000000000
5.630930319126 5.521640894956 4.450356541303 1 1 1 mag 2.000000000000 2.000000000000 2.000000000000
5.499851012568 4.003388899277 5.342621842622 1 1 1 mag 3.000000000000 3.000000000000 3.000000000000
"""
        self.assertTrue(stru_ref in c)


class TestABACUSParseStru(unittest.TestCase):
    def test_parse_stru_post(self):
        pos, move, velocity, magmom, angle1, angle2, constrain, lambda1 = (
            parse_stru_pos(
                "1.0 2.0 3.0 1 1 1 mag 1.0 2.0 3.0 v 1 1 1 angle1 100 angle2 90 sc 1 0 1 lambda 0.1 0.2 0.3"
            )
        )
        self.assertEqual(pos, [1.0, 2.0, 3.0])
        self.assertEqual(move, [1, 1, 1])
        self.assertEqual(velocity, [1.0, 1.0, 1.0])
        self.assertEqual(magmom, [1.0, 2.0, 3.0])
        self.assertEqual(angle1, 100)
        self.assertEqual(angle2, 90)
        self.assertEqual(constrain, [1, 0, 1])
        self.assertEqual(lambda1, [0.1, 0.2, 0.3])

        pos, move, velocity, magmom, angle1, angle2, constrain, lambda1 = (
            parse_stru_pos("1 2 3 mag 1 sc 1 lambda 0.1")
        )
        self.assertEqual(pos, [1, 2, 3])
        self.assertEqual(move, None)
        self.assertEqual(velocity, None)
        self.assertEqual(magmom, 1.0)
        self.assertEqual(angle1, None)
        self.assertEqual(angle2, None)
        self.assertEqual(constrain, 1)
        self.assertEqual(lambda1, 0.1)

    def test_parse_stru_error(self):
        line = "1.0 2.0 3.0 1 1"
        self.assertRaises(RuntimeError, parse_stru_pos, line), line

        line = "1.0 2.0 3.0 1 1 1 mag 1.0 3.0 v 1 1 1"
        self.assertRaises(RuntimeError, parse_stru_pos, line), line

        line = "1.0 2.0 3.0 1 1 1 mag 1 2 3 4"
        self.assertRaises(RuntimeError, parse_stru_pos, line), line

        line = "1.0 2.0 3.0 1 1 1 v 1"
        self.assertRaises(RuntimeError, parse_stru_pos, line), line

        line = "1.0 2.0 3.0 1 1 1 v 1 1"
        self.assertRaises(RuntimeError, parse_stru_pos, line), line

        line = "1.0 2.0 3.0 1 1 1 v 1 1 1 1"
        self.assertRaises(RuntimeError, parse_stru_pos, line), line

        line = "1.0 2.0 3.0 1 1 1 1"
        self.assertRaises(RuntimeError, parse_stru_pos, line), line

        line = "1.0 2.0 3.0 1 1 1 angle1 "
        self.assertRaises(RuntimeError, parse_stru_pos, line), line

        line = "1.0 2.0 3.0 1 1 1 angle1 1 2"
        self.assertRaises(RuntimeError, parse_stru_pos, line), line

        line = "1.0 2.0 3.0 1 1 1 angle2"
        self.assertRaises(RuntimeError, parse_stru_pos, line), line

        line = "1.0 2.0 3.0 angle2 1 2"
        self.assertRaises(RuntimeError, parse_stru_pos, line), line

        line = "1.0 2.0 3.0 sc"
        self.assertRaises(RuntimeError, parse_stru_pos, line), line

        line = "1.0 2.0 3.0 sc 1 2"
        self.assertRaises(RuntimeError, parse_stru_pos, line), line

        line = "1.0 2.0 3.0 lambda"
        self.assertRaises(RuntimeError, parse_stru_pos, line), line

        line = "1.0 2.0 3.0 lambda 1 2"
        self.assertRaises(RuntimeError, parse_stru_pos, line), line

        line = "1.0 2.0 3.0 lambda 1 2 3 4"
        self.assertRaises(RuntimeError, parse_stru_pos, line), line


if __name__ == "__main__":
    unittest.main()
