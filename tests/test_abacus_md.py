import os
import unittest

import numpy as np
from context import dpdata

from dpdata.unit import LengthConversion

bohr2ang = LengthConversion("bohr", "angstrom").value()


class TestABACUSMD(unittest.TestCase):
    def setUp(self):
        self.system_water = dpdata.LabeledSystem(
            "abacus.md", fmt="abacus/md"
        )  # system with stress
        self.system_Si = dpdata.LabeledSystem(
            "abacus.md.nostress", fmt="abacus/md"
        )  # system without stress
        self.system_water_unconv = dpdata.LabeledSystem(
            "abacus.md.unconv", fmt="abacus/md"
        )  # system with unconverged SCF
        self.system_newversion = dpdata.LabeledSystem(
            "abacus.md.newversion", fmt="abacus/md"
        )  # system with unconverged SCF

    def tearDown(self):
        if os.path.isfile("abacus.md/water_stru"):
            os.remove("abacus.md/water_stru")

    def test_atom_names(self):
        self.assertEqual(self.system_water.data["atom_names"], ["H", "O"])
        self.assertEqual(self.system_Si.data["atom_names"], ["Si"])
        self.assertEqual(self.system_water_unconv.data["atom_names"], ["H", "O"])

    def test_atom_numbs(self):
        self.assertEqual(self.system_water.data["atom_numbs"], [2, 1])
        self.assertEqual(self.system_Si.data["atom_numbs"], [2])
        self.assertEqual(self.system_water_unconv.data["atom_numbs"], [2, 1])

    def test_atom_types(self):
        ref_type = [0, 0, 1]
        ref_type = np.array(ref_type)
        ref_type2 = np.array([0, 0])
        ref_type3 = np.array([0, 0, 1])
        for ii in range(ref_type.shape[0]):
            self.assertEqual(self.system_water.data["atom_types"][ii], ref_type[ii])
        for ii in range(ref_type2.shape[0]):
            self.assertEqual(self.system_Si.data["atom_types"][ii], ref_type2[ii])
        for ii in range(ref_type3.shape[0]):
            self.assertEqual(
                self.system_water_unconv.data["atom_types"][ii], ref_type3[ii]
            )

    def test_cell(self):
        cell = bohr2ang * 28 * np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        cell2 = bohr2ang * 5.1 * np.array([[1, 1, 0], [1, 0, 1], [0, 1, 1]])
        cell3 = np.array(
            [
                [1.45245092e01, 0, 0],
                [-1.40550526e-02, 1.51277202e01, 0],
                [-4.42369435e-01, 4.17648184e-01, 1.49535208e01],
            ]
        )
        cell4 = np.array(
            [
                [1.24112058855e01, 0, 0],
                [0, 1.24112058855e01, 0],
                [0, 0, 1.24112058855e01],
            ]
        )
        for idx in range(np.shape(self.system_water.data["cells"])[0]):
            np.testing.assert_almost_equal(
                cell, self.system_water.data["cells"][idx], decimal=5
            )
        for idx in range(np.shape(self.system_Si.data["cells"])[0]):
            np.testing.assert_almost_equal(
                self.system_Si.data["cells"][idx], cell2, decimal=5
            )
        for idx in range(np.shape(self.system_water_unconv.data["cells"])[0]):
            np.testing.assert_almost_equal(
                self.system_water_unconv.data["cells"][idx], cell3, decimal=5
            )
        for idx in range(np.shape(self.system_newversion.data["cells"])[0]):
            np.testing.assert_almost_equal(
                self.system_newversion.data["cells"][idx], cell4, decimal=5
            )

    def test_coord(self):
        with open("abacus.md/water_coord") as fp:
            coord = []
            for ii in fp:
                coord.append([float(jj) for jj in ii.split()])
            coord = np.array(coord)
            coord = coord.reshape([5, 3, 3])
            np.testing.assert_almost_equal(
                self.system_water.data["coords"], coord, decimal=5
            )

        with open("abacus.md.nostress/Si_coord") as fp2:
            coord = []
            for ii in fp2:
                coord.append([float(jj) for jj in ii.split()])
            coord = np.array(coord)
            coord = coord.reshape([4, 2, 3])
            np.testing.assert_almost_equal(
                self.system_Si.data["coords"], coord, decimal=5
            )

        with open("abacus.md.unconv/water_coord") as fp3:
            coord = []
            for ii in fp3:
                coord.append([float(jj) for jj in ii.split()])
            coord = np.array(coord)
            coord = coord.reshape([10, 3, 3])
            np.testing.assert_almost_equal(
                self.system_water_unconv.data["coords"], coord, decimal=5
            )

        with open("abacus.md.newversion/coord.ref") as fp4:
            coord = []
            for ii in fp4:
                coord.append([float(jj) for jj in ii.split()])
            coord = np.array(coord)
            coord = coord.reshape([11, 64, 3])
            np.testing.assert_almost_equal(
                self.system_newversion.data["coords"], coord, decimal=5
            )

    def test_force(self):
        with open("abacus.md/water_force") as fp:
            force = []
            for ii in fp:
                force.append([float(jj) for jj in ii.split()])
            force = np.array(force)
            force = force.reshape([5, 3, 3])
            np.testing.assert_almost_equal(
                self.system_water.data["forces"], force, decimal=5
            )

        with open("abacus.md.nostress/Si_force") as fp2:
            force = []
            for ii in fp2:
                force.append([float(jj) for jj in ii.split()])
            force = np.array(force)
            force = force.reshape([4, 2, 3])
            np.testing.assert_almost_equal(
                self.system_Si.data["forces"], force, decimal=5
            )

        with open("abacus.md.unconv/water_force") as fp3:
            force = []
            for ii in fp3:
                force.append([float(jj) for jj in ii.split()])
            force = np.array(force)
            force = force.reshape([10, 3, 3])
            np.testing.assert_almost_equal(
                self.system_water_unconv.data["forces"], force, decimal=5
            )

        with open("abacus.md.newversion/force.ref") as fp4:
            force = []
            for ii in fp4:
                force.append([float(jj) for jj in ii.split()])
            force = np.array(force)
            force = force.reshape([11, 64, 3])
            np.testing.assert_almost_equal(
                self.system_newversion.data["forces"], force, decimal=5
            )

    def test_virial(self):
        with open("abacus.md/water_virial") as fp:
            virial = []
            for ii in fp:
                virial.append([float(jj) for jj in ii.split()])
            virial = np.array(virial)
            virial = virial.reshape([5, 3, 3])
            np.testing.assert_almost_equal(
                self.system_water.data["virials"], virial, decimal=5
            )

        with open("abacus.md.unconv/water_virial") as fp:
            virial = []
            for ii in fp:
                virial.append([float(jj) for jj in ii.split()])
            virial = np.array(virial)
            virial = virial.reshape([10, 3, 3])
            np.testing.assert_almost_equal(
                self.system_water_unconv.data["virials"], virial, decimal=5
            )

        with open("abacus.md.newversion/virial.ref") as fp:
            virial = []
            for ii in fp:
                virial.append([float(jj) for jj in ii.split()])
            virial = np.array(virial)
            virial = virial.reshape([11, 3, 3])
            np.testing.assert_almost_equal(
                self.system_newversion.data["virials"], virial, decimal=5
            )

    def test_energy(self):
        ref_energy = np.array(
            [-466.69285117, -466.69929051, -466.69829826, -466.70364664, -466.6976083]
        )
        ref_energy2 = np.array(
            [-211.77184603, -211.78111966, -211.79681663, -211.79875524]
        )
        ref_energy3 = np.array(
            [
                -464.87380991,
                -465.18489358,
                -465.97407849,
                -465.98292836,
                -465.85528926,
                -465.33957501,
                -464.64886682,
                -464.61802032,
                -465.61854656,
                -466.05660096,
            ]
        )
        np.testing.assert_almost_equal(self.system_water.data["energies"], ref_energy)
        np.testing.assert_almost_equal(self.system_Si.data["energies"], ref_energy2)
        np.testing.assert_almost_equal(
            self.system_water_unconv.data["energies"], ref_energy3
        )

    def test_to_system(self):
        pp_file = ["H.upf", "O.upf"]
        numerical_orbital = ["H.upf", "O.upf"]
        numerical_descriptor = "jle.orb"
        mass = [1.008, 15.994]
        self.system_water.to(
            file_name="abacus.md/water_stru",
            fmt="abacus/stru",
            pp_file=pp_file,
            numerical_orbital=numerical_orbital,
            numerical_descriptor=numerical_descriptor,
            mass=mass,
        )
        self.assertTrue(os.path.isfile("abacus.md/water_stru"))
        if os.path.isfile("abacus.md/water_stru"):
            with open("abacus.md/water_stru") as f:
                iline = 0
                for iline, l in enumerate(f):
                    iline += 1
                self.assertEqual(iline, 30)


if __name__ == "__main__":
    unittest.main()
