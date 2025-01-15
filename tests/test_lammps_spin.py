from __future__ import annotations

import os
import shutil
import unittest

import numpy as np
from context import dpdata

from dpdata.lammps.dump import get_spin

TRAJ_NO_ID = """ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
17
ITEM: BOX BOUNDS xy xz yz pp pp pp
-4.0080511965879438e-02 5.7039029418994556e+00 -5.9179115295410201e-03
1.4436085788922526e-02 5.6674744441011660e+00 -1.1487414836883500e-02
7.8239288740356017e-03 5.6734038274259646e+00 6.8277359008788905e-04
ITEM: ATOMS type x y z c_spin[1] c_spin[2] c_spin[3] c_spin[4] c_spin[5] c_spin[6] c_spin[7] c_spin[8] c_spin[9] c_spin[10]
1 0.00141160 5.64868599 0.01005602 1.54706291 0.00000000 0.00000000 1.00000000 -1.40772100 -2.03739417 -1522.64797384 -0.00397809 -0.00190426 -0.00743976
1 5.65283939 5.57449025 2.84281508 1.54412869 0.00000000 0.00000000 1.00000000 7.75304092 6.48949619 -1512.84926162 -0.00637234 -0.00733168 0.00661107
1 0.00066480 2.78022036 0.01010716 1.54612979 0.00000000 0.00000000 1.00000000 -0.93618575 1.92206111 -1520.80305011 -0.00316673 0.00177893 -0.00744575
1 5.65233666 2.85374747 2.84289453 1.54439093 0.00000000 0.00000000 1.00000000 8.11012818 -6.49922039 -1514.31557088 -0.00569217 0.00741000 0.00640353
1 2.82063515 5.64869321 0.01007552 1.54714250 0.00000000 0.00000000 1.00000000 2.49070852 -2.14456666 -1523.53038650 0.00478410 -0.00213962 -0.00751154
1 2.89579803 5.57439179 2.84287630 1.54415032 0.00000000 0.00000000 1.00000000 -8.03062338 6.63950296 -1513.41291897 0.00440396 -0.00717185 0.00633657
1 2.82151287 2.78010538 0.01016303 1.54619615 0.00000000 0.00000000 1.00000000 2.71859584 1.98482729 -1521.34149633 0.00533453 0.00194532 -0.00745901
1 2.89637049 2.85377083 2.84297332 1.54440023 0.00000000 0.00000000 1.00000000 -7.76758760 -6.67134514 -1514.43304618 0.00505040 0.00743195 0.00630302
1 1.41106492 1.38817482 1.72302072 1.18134529 0.00000000 0.00000000 1.00000000 0.27170165 -0.00426695 -444.22843899 0.00100237 -0.00002725 -0.03385965
 1 1.41105247 1.38807861 3.96314606 1.18153407 0.00000000 0.00000000 1.00000000 -0.07722674 0.01368756 -337.08703133 -0.00066982 0.00007487 0.07887183
 1 1.41105864 4.21395432 1.43987180 1.71989299 0.00000000 0.00000000 1.00000000 -0.01511106 0.00320081 -1653.34500916 0.00010421 0.00007248 0.00634401
 1 1.41104843 4.21387554 4.24576823 1.71989825 0.00000000 0.00000000 1.00000000 -0.71645898 0.05923960 -1640.68070568 -0.00117959 0.00006676 -0.01467806
 1 4.27433865 1.38779084 1.43977211 1.72010048 0.00000000 0.00000000 1.00000000 0.45899480 0.03956420 -1653.36356942 0.00051885 0.00002313 0.00911600
 1 4.27436799 1.38772964 4.24586490 1.72010133 0.00000000 0.00000000 1.00000000 0.38385331 0.07301994 -1642.06086017 -0.00002034 0.00010335 -0.01688908
 1 4.27435427 4.21452597 1.39359689 1.65590121 0.00000000 0.00000000 1.00000000 -0.01658773 -0.06159007 -1659.12744163 0.00006470 -0.00006420 -0.01342935
 1 4.27434583 4.21455469 4.29208004 1.65592002 0.00000000 0.00000000 1.00000000 -0.15590720 -0.03252166 -1654.84697132 -0.00066755 -0.00003915 -0.00482188
 2 1.41096761 1.38958048 0.01029027 0.00000000 0.00000000 0.00000000 1.00000000 0.00000000 0.00000000 0.00000000 0.00048351 -0.00022876 -0.00645195"""


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

    def test_read_dump_partial_spin(self):
        # test if dpdata can read the spins when the spin data is not complete
        with self.assertWarns(UserWarning) as cm:
            tmp_system = dpdata.System(
                "lammps/traj_partial_spin.dump",
                fmt="lammps/dump",
                type_map=["O", "H"],
                input_file="lammps/in.lmp",
            )
            self.assertTrue("spins" not in tmp_system.data)

        self.assertIn("Warning: spin info is not found in frame", str(cm.warning))

    def test_get_spin_failed(self):
        with self.assertWarns(UserWarning) as cm:
            spin = get_spin(
                TRAJ_NO_ID.split("\n"),
                ["c_spin[1]", "c_spin[2]", "c_spin[3]", "c_spin[4]"],
            )
            self.assertTrue(spin is None)

        self.assertIn("Error processing spin data:", str(cm.warning))
