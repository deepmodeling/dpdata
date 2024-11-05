from __future__ import annotations

import os
import shutil
import unittest

import numpy as np
from context import dpdata


class TestABACUSSpin(unittest.TestCase):
    def setUp(self):
        self.dump_path = "abacus.spin/dump"
        os.makedirs(self.dump_path, exist_ok=True)

    def tearDown(self):
        if os.path.isdir(self.dump_path):
            shutil.rmtree(self.dump_path)

    def test_scf(self):
        os.system("cp abacus.spin/INPUT.scf abacus.spin/INPUT")
        mysys = dpdata.LabeledSystem("abacus.spin", fmt="abacus/scf")
        data = mysys.data
        self.assertAlmostEqual(data["energies"][0], -6818.719409466637)
        np.testing.assert_almost_equal(
            data["spins"][0],
            [
                [-0.0000002724, -0.0000001728, 2.4000001004],
                [-0.0000003180, -0.0000002299, 2.3999994597],
            ],
            decimal=8,
        )
        np.testing.assert_almost_equal(
            data["force_mags"][0],
            [
                [-0.0000175013, -0.0000418680, -0.3669618965],
                [-0.0000161517, -0.0000195198, -0.3669821632],
            ],
            decimal=8,
        )

        # dump to deepmd-npy
        mysys.to(file_name=self.dump_path, fmt="deepmd/npy")
        self.assertTrue(os.path.isfile(f"{self.dump_path}/set.000/spin.npy"))
        self.assertTrue(os.path.isfile(f"{self.dump_path}/set.000/force_mag.npy"))

        sys2 = dpdata.LabeledSystem(self.dump_path, fmt="deepmd/npy")
        np.testing.assert_almost_equal(data["spins"], sys2.data["spins"], decimal=8)
        np.testing.assert_almost_equal(
            data["force_mags"], sys2.data["force_mags"], decimal=8
        )

    def test_relax(self):
        os.system("cp abacus.spin/INPUT.relax abacus.spin/INPUT")
        mysys = dpdata.LabeledSystem("abacus.spin", fmt="abacus/relax")
        data = mysys.data
        spins_ref = np.array(
            [
                [
                    [1.16909819, 1.16895965, 1.16895485],
                    [1.16827825, 1.16832716, 1.16836899],
                ],
                [
                    [1.25007143, 1.25006167, 1.25004587],
                    [1.25015764, 1.2501678, 1.25018344],
                ],
                [
                    [1.24984994, 1.24977108, 1.24978313],
                    [1.24996533, 1.2500441, 1.25003208],
                ],
            ]
        )
        magforces_ref = np.array(
            [
                [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                [
                    [-0.16734626, -0.16735378, -0.16735617],
                    [-0.16836467, -0.16835897, -0.16835625],
                ],
                [
                    [-0.16573406, -0.16574627, -0.1657445],
                    [-0.16619489, -0.16617948, -0.16618272],
                ],
            ]
        )
        self.assertEqual(len(data["spins"]), 3)
        self.assertEqual(len(data["force_mags"]), 3)
        np.testing.assert_almost_equal(data["spins"], spins_ref, decimal=8)
        np.testing.assert_almost_equal(data["force_mags"], magforces_ref, decimal=8)

        # dump to deepmd-npy
        mysys.to(file_name=self.dump_path, fmt="deepmd/npy")
        self.assertTrue(os.path.isfile(f"{self.dump_path}/set.000/spin.npy"))
        self.assertTrue(os.path.isfile(f"{self.dump_path}/set.000/force_mag.npy"))

        sys2 = dpdata.LabeledSystem(self.dump_path, fmt="deepmd/npy")
        np.testing.assert_almost_equal(data["spins"], sys2.data["spins"], decimal=8)
        np.testing.assert_almost_equal(
            data["force_mags"], sys2.data["force_mags"], decimal=8
        )

    def test_md(self):
        os.system("cp abacus.spin/INPUT.md abacus.spin/INPUT")
        mysys = dpdata.LabeledSystem("abacus.spin", fmt="abacus/md")
        data = mysys.data
        spins_ref = np.array(
            [
                [
                    [1.16909819, 1.16895965, 1.16895485],
                    [1.16827825, 1.16832716, 1.16836899],
                ],
                [
                    [1.2500362, 1.25007501, 1.2500655],
                    [1.25019078, 1.25015253, 1.25016188],
                ],
                [
                    [1.24985138, 1.24976901, 1.2497695],
                    [1.24996388, 1.25004618, 1.25004561],
                ],
                [
                    [1.24982513, 1.24985445, 1.24985336],
                    [1.25005073, 1.25001814, 1.25002065],
                ],
            ]
        )
        magforces_ref = np.array(
            [
                [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                [
                    [-0.16747275, -0.16747145, -0.16746776],
                    [-0.16853881, -0.16853935, -0.16854119],
                ],
                [
                    [-0.16521817, -0.16523256, -0.16523212],
                    [-0.16549418, -0.16547867, -0.16547913],
                ],
                [
                    [-0.16141172, -0.16140644, -0.1614127],
                    [-0.15901519, -0.15905932, -0.15904824],
                ],
            ]
        )
        self.assertEqual(len(data["spins"]), 4)
        self.assertEqual(len(data["force_mags"]), 4)
        np.testing.assert_almost_equal(data["spins"], spins_ref, decimal=8)
        np.testing.assert_almost_equal(data["force_mags"], magforces_ref, decimal=8)

        # dump to deepmd-npy
        mysys.to(file_name=self.dump_path, fmt="deepmd/npy")
        self.assertTrue(os.path.isfile(f"{self.dump_path}/set.000/spin.npy"))
        self.assertTrue(os.path.isfile(f"{self.dump_path}/set.000/force_mag.npy"))

        sys2 = dpdata.LabeledSystem(self.dump_path, fmt="deepmd/npy")
        np.testing.assert_almost_equal(data["spins"], sys2.data["spins"], decimal=8)
        np.testing.assert_almost_equal(
            data["force_mags"], sys2.data["force_mags"], decimal=8
        )
