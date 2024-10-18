from __future__ import annotations

import os
import shutil
import unittest

import numpy as np
from context import dpdata


class TestLmp(unittest.TestCase):
    def setUp(self):
        self.tmp_system = dpdata.System(
            os.path.join("poscars", "conf.lmp"), type_map=["O", "H"]
        )
        self.tmp_system.data["spins"] = [[[3, 4, 0], [0, 4, 3]]]
        self.lmp_coord_name = "tmp.lmp"

    def tearDown(self):
        pass  # if os.path.isfile(self.lmp_coord_name):os.remove(self.lmp_coord_name)

    def test_dump_input(self):
        self.tmp_system.to("lammps/lmp", self.lmp_coord_name)
        self.assertTrue(os.path.isfile(self.lmp_coord_name))
        with open(self.lmp_coord_name) as f:
            c = f.read()

        coord_ref = """     1      1    0.0000000000    0.0000000000    0.0000000000    0.6000000000    0.8000000000    0.0000000000    5.0000000000
     2      2    1.2621856000    0.7018028000    0.5513885000    0.0000000000    0.8000000000    0.6000000000    5.0000000000"""
        self.assertTrue(coord_ref in c)

    def test_dump_input_zero_spin(self):
        self.tmp_system.data["spins"] = [[[0, 0, 0], [0, 0, 0]]]
        self.tmp_system.to("lammps/lmp", self.lmp_coord_name)
        self.assertTrue(os.path.isfile(self.lmp_coord_name))
        with open(self.lmp_coord_name) as f:
            c = f.read()
        coord_ref = """     1      1    0.0000000000    0.0000000000    0.0000000000    0.0000000000    0.0000000000    1.0000000000    0.0000000000
     2      2    1.2621856000    0.7018028000    0.5513885000    0.0000000000    0.0000000000    1.0000000000    0.0000000000"""
        self.assertTrue(coord_ref in c)

    def test_read_input(self):
        # check if dpdata can read the spins
        tmp_system = dpdata.System(
            "lammps/spin.lmp", fmt="lammps/lmp", type_map=["O", "H"]
        )
        self.assertTrue((tmp_system.data["spins"][0] == [[3, 4, 0], [0, 4, 3]]).all())

        tmp_system.to(file_name="lammps/dump", fmt="deepmd/npy")
        self.assertTrue(os.path.isfile("lammps/dump/set.000/spin.npy"))

        if os.path.isdir("lammps/dump"):
            shutil.rmtree("lammps/dump")


class TestDump(unittest.TestCase):
    def test_read_dump_spin(self):
        tmp_system = dpdata.System(
            "lammps/traj.dump",
            fmt="lammps/dump",
            type_map=["O", "H"],
            input_file="lammps/in.lmp",
        )
        self.assertTrue(len(tmp_system.data["spins"]) == 2)
        np.testing.assert_almost_equal(
            tmp_system.data["spins"][0][0], [0, 0, 1.54706291], decimal=8
        )
        np.testing.assert_almost_equal(
            tmp_system.data["spins"][0][1], [0, 0, 1.54412869], decimal=8
        )
        np.testing.assert_almost_equal(
            tmp_system.data["spins"][0][-2], [0, 0, 1.65592002], decimal=8
        )
        np.testing.assert_almost_equal(
            tmp_system.data["spins"][0][-1], [0, 0, 0], decimal=8
        )

        np.testing.assert_almost_equal(
            tmp_system.data["spins"][1][0],
            [0.21021514724299958, 1.0123821159859323, -0.6159960941686954],
            decimal=8,
        )
        np.testing.assert_almost_equal(
            tmp_system.data["spins"][1][1],
            [1.0057302798645609, 0.568273899191638, -0.2363447073875224],
            decimal=8,
        )
        np.testing.assert_almost_equal(
            tmp_system.data["spins"][1][-2],
            [-0.28075943761984146, -1.2845200151690905, -0.0201237855118935],
            decimal=8,
        )
        np.testing.assert_almost_equal(
            tmp_system.data["spins"][1][-1], [0, 0, 0], decimal=8
        )

        tmp_system.to(file_name="lammps/dump", fmt="deepmd/npy")
        self.assertTrue(os.path.isfile("lammps/dump/set.000/spin.npy"))

        if os.path.isdir("lammps/dump"):
            shutil.rmtree("lammps/dump")
