from __future__ import annotations

import os
import unittest

import numpy as np
from comp_sys import CompSys, IsPBC
from context import dpdata
from poscars.poscar_ref_oh import TestPOSCARoh


class TestPOSCARCart(unittest.TestCase, TestPOSCARoh):
    def setUp(self):
        self.system = dpdata.System()
        self.system.from_vasp_poscar(os.path.join("poscars", "POSCAR.oh.c"))

    def test_move_flags(self):
        expected = np.array([[[True, True, False], [False, False, False]]])
        self.assertTrue(np.array_equal(self.system["move"], expected))


class TestPOSCARMoveFlags(unittest.TestCase):
    def test_move_flags_error1(self):
        with self.assertRaisesRegex(RuntimeError, "Invalid move flags.*?"):
            dpdata.System().from_vasp_poscar(os.path.join("poscars", "POSCAR.oh.err1"))

    def test_move_flags_error2(self):
        with self.assertRaisesRegex(RuntimeError, "Invalid move flag: a"):
            dpdata.System().from_vasp_poscar(os.path.join("poscars", "POSCAR.oh.err2"))

    def test_move_flags_error3(self):
        system = dpdata.System().from_vasp_poscar(
            os.path.join("poscars", "POSCAR.oh.c")
        )
        system.data["move"] = np.array([[[True, True], [False, False]]])
        with self.assertRaisesRegex(
            RuntimeError, "Invalid move flags:.*?should be a list of 3 bools"
        ):
            system.to_vasp_poscar("POSCAR.tmp.1")


class TestPOSCARDirect(unittest.TestCase, TestPOSCARoh):
    def setUp(self):
        self.system = dpdata.System()
        self.system.from_vasp_poscar(os.path.join("poscars", "POSCAR.oh.d"))


class TestPOSCARDirectDuplicated(unittest.TestCase):
    def test(self):
        ss = dpdata.System(
            os.path.join("poscars", "POSCAR.oh.d.dup"), fmt="vasp/poscar"
        )
        self.assertEqual(ss["atom_names"], ["O", "H"])
        self.assertEqual(ss["atom_numbs"], [2, 1])
        self.assertEqual(list(ss["atom_types"]), [0, 1, 0])

    def test_type_map(self):
        ss = dpdata.System(
            os.path.join("poscars", "POSCAR.oh.d.dup"),
            fmt="vasp/poscar",
            type_map=["H", "O"],
        )
        self.assertEqual(ss["atom_names"], ["H", "O"])
        self.assertEqual(ss["atom_numbs"], [1, 2])
        self.assertEqual(list(ss["atom_types"]), [1, 0, 1])


class TestVaspPOSCARTypeMap(unittest.TestCase, CompSys, IsPBC):
    def setUp(self):
        sys0 = dpdata.System("poscars/POSCAR.oh.d", fmt="vasp/poscar")
        sys0.data["atom_names"] = ["A", "H", "B", "O", "D"]
        sys0.data["atom_numbs"] = [0, 1, 0, 1, 0]
        sys0.data["atom_types"] = np.array([3, 1], dtype=int)
        sys1 = dpdata.System(
            "poscars/POSCAR.oh.d", fmt="vasp/poscar", type_map=["A", "H", "B", "O", "D"]
        )
        self.system_1 = sys0
        self.system_2 = sys1
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6


if __name__ == "__main__":
    unittest.main()
