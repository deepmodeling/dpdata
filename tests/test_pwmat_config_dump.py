from __future__ import annotations

import os
import unittest

import numpy as np
from pwmat.config_ref_oh import Testconfigoh

import dpdata
from dpdata.formats.pwmat.atomconfig import from_system_data


def myfilecmp(test, f0, f1):
    with open(f0) as fp0:
        with open(f1) as fp1:
            test.assertTrue(fp0.read() == fp1.read())


class TestatomconfigDump(unittest.TestCase, Testconfigoh):
    def setUp(self):
        tmp_system = dpdata.System()
        tmp_system.from_lammps_lmp(
            os.path.join("pwmat", "conf.lmp"), type_map=["O", "H"]
        )
        tmp_system.to_pwmat_atomconfig("tmp.atom.config")
        self.system = dpdata.System()
        self.system.from_pwmat_atomconfig("tmp.atom.config")


class TestatomconfigDump1(unittest.TestCase, Testconfigoh):
    def setUp(self):
        tmp_system = dpdata.System()
        tmp_system.from_pwmat_atomconfig(os.path.join("pwmat", "atom.config.oh"))
        # tmp_system.from_lammps_lmp(os.path.join('poscars', 'conf.lmp'), type_map = ['O', 'H'])
        tmp_system.to_pwmat_atomconfig("tmp.atom.config")
        self.system = dpdata.System()
        self.system.from_pwmat_atomconfig("tmp.atom.config")


class TestatomconfigSkipZeroAtomNumb(unittest.TestCase):
    def tearDown(self):
        if os.path.isfile("atom.config.tmp.1"):
            os.remove("atom.config.tmp.1")
        if os.path.isfile("atom.config.tmp.2"):
            os.remove("atom.config.tmp.2")

    def test_dump_pwmat_type_map(self):
        system0 = dpdata.System(
            os.path.join("pwmat", "atom.config.oh"),
            fmt="pwmat/atom.config",
            type_map=["H", "O"],
        )
        system0.to_pwmat_atomconfig("atom.config.tmp.1")
        system1 = dpdata.System(
            os.path.join("pwmat", "atom.config.oh"),
            fmt="pwmat/atom.config",
            type_map=["C", "H", "A", "O", "B"],
        )
        system1.to_pwmat_atomconfig("atom.config.tmp.2")
        myfilecmp(self, "atom.config.tmp.1", "atom.config.tmp.2")


class TestAtomconfigFrameCell(unittest.TestCase):
    def test_fractional_coordinates_use_selected_frame_cell(self):
        system = {
            "atom_names": ["H"],
            "atom_numbs": [1],
            "atom_types": np.array([0]),
            "cells": np.array([np.eye(3), np.eye(3) * 2.0]),
            "coords": np.array([[[1.0, 0.0, 0.0]], [[1.0, 0.0, 0.0]]]),
        }
        output = from_system_data(system, f_idx=1)
        self.assertIn("0.5000000000    0.0000000000    0.0000000000", output)


if __name__ == "__main__":
    unittest.main()
