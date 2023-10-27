import os
import shutil
import unittest

import numpy as np
from context import dpdata

from dpdata.unit import LengthConversion

bohr2ang = LengthConversion("bohr", "angstrom").value()


class TestABACUSRelaxLabeledOutput(unittest.TestCase):
    def setUp(self):
        shutil.copy(
            "abacus.relax/OUT.abacus/running_cell-relax.log.normal",
            "abacus.relax/OUT.abacus/running_cell-relax.log",
        )
        self.system = dpdata.LabeledSystem("abacus.relax", fmt="abacus/relax")

    def tearDown(self):
        if os.path.isfile("abacus.relax/OUT.abacus/running_cell-relax.log"):
            os.remove("abacus.relax/OUT.abacus/running_cell-relax.log")

    def test_atom_names(self):
        self.assertEqual(self.system.data["atom_names"], ["H", "O"])

    def test_atom_numbs(self):
        self.assertEqual(self.system.data["atom_numbs"], [2, 1])

    def test_atom_types(self):
        ref_type = np.array([0, 0, 1])
        for ii in range(ref_type.shape[0]):
            self.assertEqual(self.system.data["atom_types"][ii], ref_type[ii])

    def test_cell(self):
        cell = bohr2ang * 28.0 * np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        for idx in range(np.shape(self.system.data["cells"])[0]):
            np.testing.assert_almost_equal(
                cell, self.system.data["cells"][idx], decimal=5
            )

    def test_coord(self):
        with open("abacus.relax/coord.ref") as fp:
            ref = []
            for ii in fp:
                ref.append([float(jj) for jj in ii.split()])
            ref = np.array(ref)
            ref = ref.reshape([5, 3, 3])
            np.testing.assert_almost_equal(self.system.data["coords"], ref, decimal=5)

    def test_force(self):
        with open("abacus.relax/force.ref") as fp:
            ref = []
            for ii in fp:
                ref.append([float(jj) for jj in ii.split()])
            ref = np.array(ref)
            ref = ref.reshape([5, 3, 3])
            np.testing.assert_almost_equal(self.system.data["forces"], ref, decimal=5)

    def test_virial(self):
        with open("abacus.relax/virial.ref") as fp:
            ref = []
            for ii in fp:
                ref.append([float(jj) for jj in ii.split()])
            ref = np.array(ref)
            ref = ref.reshape([5, 3, 3])
            np.testing.assert_almost_equal(self.system.data["virials"], ref, decimal=5)

    def test_stress(self):
        with open("abacus.relax/stress.ref") as fp:
            ref = []
            for ii in fp:
                ref.append([float(jj) for jj in ii.split()])
            ref = np.array(ref)
            ref = ref.reshape([5, 3, 3])
            np.testing.assert_almost_equal(self.system.data["stress"], ref, decimal=5)

    def test_energy(self):
        ref_energy = np.array(
            [-465.77753104, -464.35757552, -465.79307346, -465.80056811, -465.81235433]
        )
        np.testing.assert_almost_equal(self.system.data["energies"], ref_energy)


class TestABACUSRelaxLabeledOutputAbnormal(unittest.TestCase):
    def setUp(self):
        shutil.copy(
            "abacus.relax/OUT.abacus/running_cell-relax.log.abnormal",
            "abacus.relax/OUT.abacus/running_cell-relax.log",
        )
        self.system = dpdata.LabeledSystem("abacus.relax", fmt="abacus/relax")

    def test_result(self):
        data = self.system.data
        self.assertEqual(len(data["coords"]), 4)
        self.assertEqual(len(data["energies"]), len(data["coords"]))
        self.assertEqual(len(data["cells"]), len(data["coords"]))
        self.assertEqual(len(data["forces"]), len(data["coords"]))
        self.assertEqual(len(data["stress"]), len(data["coords"]))
        self.assertEqual(len(data["virials"]), len(data["coords"]))
        np.testing.assert_almost_equal(data["energies"][3], -465.81235433)

    def tearDown(self):
        if os.path.isfile("abacus.relax/OUT.abacus/running_cell-relax.log"):
            os.remove("abacus.relax/OUT.abacus/running_cell-relax.log")


class TestABACUSRelaxLabeledOutputNoStress(unittest.TestCase):
    def setUp(self):
        shutil.copy(
            "abacus.relax/OUT.abacus/running_cell-relax.log.nostress",
            "abacus.relax/OUT.abacus/running_cell-relax.log",
        )
        shutil.move(
            "abacus.relax/STRU",
            "abacus.relax/STRU.bak",
        )
        shutil.copy(
            "abacus.relax/STRU.Si",
            "abacus.relax/STRU",
        )

    def test_result(self):
        system = dpdata.LabeledSystem("abacus.relax", fmt="abacus/relax")
        self.assertRaises(KeyError, lambda: system.data["virials"])

    def tearDown(self):
        if os.path.isfile("abacus.relax/OUT.abacus/running_cell-relax.log"):
            os.remove("abacus.relax/OUT.abacus/running_cell-relax.log")
        shutil.move(
            "abacus.relax/STRU.bak",
            "abacus.relax/STRU",
        )


class TestABACUSRelaxLabeledOutputV341(unittest.TestCase):
    # Since ABACUS v3.4.1, the output format of force and stress has been changed.
    def setUp(self):
        shutil.copy(
            "abacus.relax/OUT.abacus/running_cell-relax.log.v341",
            "abacus.relax/OUT.abacus/running_cell-relax.log",
        )
        shutil.move(
            "abacus.relax/STRU",
            "abacus.relax/STRU.bak",
        )
        shutil.copy(
            "abacus.relax/STRU.Si",
            "abacus.relax/STRU",
        )
        self.system = dpdata.LabeledSystem("abacus.relax", fmt="abacus/relax")

    def test_force(self):
        # with open("abacus.relax/force.v341.ref","w") as fp:
        #    for i in self.system.data["forces"]:
        #        for j in i:
        #            fp.write("%f %f %f\n" % tuple(j.tolist()))
        with open("abacus.relax/force.v341.ref") as fp:
            ref = []
            for ii in fp:
                ref.append([float(jj) for jj in ii.split()])
            ref = np.array(ref)
            ref = ref.reshape([4, 2, 3])
            np.testing.assert_almost_equal(self.system.data["forces"], ref, decimal=5)

    def test_stress(self):
        # with open("abacus.relax/stress.v341.ref","w") as fp:
        #    for i in self.system.data["stress"]:
        #        for j in i:
        #            fp.write("%f %f %f\n" % tuple(j.tolist()))
        with open("abacus.relax/stress.v341.ref") as fp:
            ref = []
            for ii in fp:
                ref.append([float(jj) for jj in ii.split()])
            ref = np.array(ref)
            ref = ref.reshape([4, 3, 3])
            np.testing.assert_almost_equal(self.system.data["stress"], ref, decimal=5)

    def tearDown(self):
        if os.path.isfile("abacus.relax/OUT.abacus/running_cell-relax.log"):
            os.remove("abacus.relax/OUT.abacus/running_cell-relax.log")
        shutil.move(
            "abacus.relax/STRU.bak",
            "abacus.relax/STRU",
        )


if __name__ == "__main__":
    unittest.main()
