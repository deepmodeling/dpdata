from __future__ import annotations

import unittest

import numpy as np
from context import dpdata


class TestPWSCFSinglePointEnergy:
    def test_atom_names(self):
        self.assertEqual(self.system_ch4.data["atom_names"], ["H", "C"])
        self.assertEqual(self.system_h2o.data["atom_names"], ["O", "H"])
        self.assertEqual(self.system_ch4_2.data["atom_names"], ["H", "C"])

    def test_atom_numbs(self):
        self.assertEqual(self.system_ch4.data["atom_numbs"], [4, 1])
        self.assertEqual(self.system_h2o.data["atom_numbs"], [64, 128])
        self.assertEqual(self.system_ch4_2.data["atom_numbs"], [4, 1])

    def test_atom_types(self):
        ref_type = [0, 0, 0, 0, 1]
        ref_type = np.array(ref_type)
        for ii in range(ref_type.shape[0]):
            self.assertEqual(self.system_ch4.data["atom_types"][ii], ref_type[ii])

        ref_type = [0] * 64 + [1] * 128
        ref_type = np.array(ref_type)
        for ii in range(ref_type.shape[0]):
            self.assertEqual(self.system_h2o.data["atom_types"][ii], ref_type[ii])

        ref_type = [0, 0, 0, 0, 1]
        ref_type = np.array(ref_type)
        for ii in range(ref_type.shape[0]):
            self.assertEqual(self.system_ch4_2.data["atom_types"][ii], ref_type[ii])

    def test_cell(self):
        cell = 10 * np.eye(3)
        for ii in range(cell.shape[0]):
            for jj in range(cell.shape[1]):
                self.assertAlmostEqual(
                    self.system_ch4.data["cells"][0][ii][jj], cell[ii][jj]
                )

        fp = open("qe.scf/h2o_cell")
        cell = []
        for ii in fp:
            cell.append([float(jj) for jj in ii.split()])
        cell = np.array(cell)
        for ii in range(cell.shape[0]):
            for jj in range(cell.shape[1]):
                self.assertAlmostEqual(
                    self.system_h2o.data["cells"][0][ii][jj], cell[ii][jj]
                )
        fp.close()

        cell = 10 * np.eye(3)
        for ii in range(cell.shape[0]):
            for jj in range(cell.shape[1]):
                self.assertAlmostEqual(
                    self.system_ch4_2.data["cells"][0][ii][jj], cell[ii][jj]
                )

    def test_coord(self):
        fp = open("qe.scf/ch4_coord")
        coord = []
        for ii in fp:
            coord.append([float(jj) for jj in ii.split()])
        coord = np.array(coord)
        for ii in range(coord.shape[0]):
            for jj in range(coord.shape[1]):
                self.assertAlmostEqual(
                    self.system_ch4.data["coords"][0][ii][jj], coord[ii][jj]
                )
        fp.close()

        fp = open("qe.scf/h2o_coord")
        coord = []
        for ii in fp:
            coord.append([float(jj) for jj in ii.split()])
        coord = np.array(coord)
        for ii in range(coord.shape[0]):
            for jj in range(coord.shape[1]):
                self.assertAlmostEqual(
                    self.system_h2o.data["coords"][0][ii][jj], coord[ii][jj]
                )
        fp.close()

        fp = open("qe.scf/ch4_coord")
        coord = []
        for ii in fp:
            coord.append([float(jj) for jj in ii.split()])
        coord = np.array(coord)
        for ii in range(coord.shape[0]):
            for jj in range(coord.shape[1]):
                self.assertAlmostEqual(
                    self.system_ch4_2.data["coords"][0][ii][jj], coord[ii][jj]
                )
        fp.close()

    def test_force(self):
        fp = open("qe.scf/ch4_force")
        force = []
        for ii in fp:
            force.append([float(jj) for jj in ii.split()])
        force = np.array(force)
        for ii in range(force.shape[0]):
            for jj in range(force.shape[1]):
                self.assertAlmostEqual(
                    self.system_ch4.data["forces"][0][ii][jj], force[ii][jj]
                )
        fp.close()

        fp = open("qe.scf/h2o_force")
        force = []
        for ii in fp:
            force.append([float(jj) for jj in ii.split()])
        force = np.array(force)
        for ii in range(force.shape[0]):
            for jj in range(force.shape[1]):
                self.assertAlmostEqual(
                    self.system_h2o.data["forces"][0][ii][jj], force[ii][jj]
                )
        fp.close()

        fp = open("qe.scf/ch4_force_2")
        force = []
        for ii in fp:
            force.append([float(jj) for jj in ii.split()])
        force = np.array(force)
        for ii in range(force.shape[0]):
            for jj in range(force.shape[1]):
                self.assertAlmostEqual(
                    self.system_ch4_2.data["forces"][0][ii][jj], force[ii][jj]
                )
        fp.close()

    def test_virial(self):
        fp = open("qe.scf/ch4_virial")
        virial = []
        for ii in fp:
            virial.append([float(jj) for jj in ii.split()])
        virial = np.array(virial)
        for ii in range(virial.shape[0]):
            for jj in range(virial.shape[1]):
                self.assertAlmostEqual(
                    self.system_ch4.data["virials"][0][ii][jj], virial[ii][jj], places=3
                )
        fp.close()

        fp = open("qe.scf/h2o_virial")
        virial = []
        for ii in fp:
            virial.append([float(jj) for jj in ii.split()])
        virial = np.array(virial)
        for ii in range(virial.shape[0]):
            for jj in range(virial.shape[1]):
                self.assertAlmostEqual(
                    self.system_h2o.data["virials"][0][ii][jj], virial[ii][jj], places=2
                )
        fp.close()

    def test_energy(self):
        ref_energy = -219.74425946528794
        self.assertAlmostEqual(self.system_ch4.data["energies"][0], ref_energy)
        ref_energy = -30007.651851226798
        self.assertAlmostEqual(self.system_h2o.data["energies"][0], ref_energy)
        ref_energy = -219.7153691367526562
        self.assertAlmostEqual(self.system_ch4_2.data["energies"][0], ref_energy)


class TestPWSCFLabeledOutput(unittest.TestCase, TestPWSCFSinglePointEnergy):
    def setUp(self):
        self.system_ch4 = dpdata.LabeledSystem("qe.scf/01.out", fmt="qe/pw/scf")
        self.system_h2o = dpdata.LabeledSystem("qe.scf/02.out", fmt="qe/pw/scf")
        self.system_ch4_2 = dpdata.LabeledSystem("qe.scf/03.out", fmt="qe/pw/scf")


class TestPWSCFLabeledOutputListInput(unittest.TestCase, TestPWSCFSinglePointEnergy):
    def setUp(self):
        self.system_ch4 = dpdata.LabeledSystem(
            ["qe.scf/01.in", "qe.scf/01.out"], fmt="qe/pw/scf"
        )
        self.system_h2o = dpdata.LabeledSystem(
            ["qe.scf/02.in", "qe.scf/02.out"], fmt="qe/pw/scf"
        )
        self.system_ch4_2 = dpdata.LabeledSystem(
            ["qe.scf/03.in", "qe.scf/03.out"], fmt="qe/pw/scf"
        )


class TestNa(unittest.TestCase):
    def test(self):
        ss = dpdata.LabeledSystem("qe.scf/na.out", fmt="qe/pw/scf")
        self.assertEqual(ss["atom_numbs"], [3])
        self.assertEqual(ss["atom_names"], ["Na"])
        self.assertEqual(ss.get_nframes(), 1)
        np.testing.assert_array_equal(ss["atom_types"], [0, 0, 0])
        np.testing.assert_allclose(
            ss["coords"][0],
            np.array(
                [
                    0.940587444301534,
                    0.397635863676890,
                    0.059472381901808,
                    0.059413515648061,
                    0.602364552456546,
                    0.559472465518034,
                    0.602364619812068,
                    0.059413062489401,
                    0.059472381901808,
                ]
            ).reshape(3, 3)
            @ ss["cells"][0],
        )
        np.testing.assert_allclose(
            ss["cells"][0],
            np.array(
                [
                    7.171683039200000,
                    0.000000000000000,
                    0.000000000000000,
                    -4.270578118300000,
                    5.761527588200000,
                    0.000000000000000,
                    -0.000000045600000,
                    0.000000023000000,
                    12.826457854999999,
                ]
            ).reshape(3, 3),
        )
        np.testing.assert_allclose(ss["forces"][0], np.zeros([3, 3]))


if __name__ == "__main__":
    unittest.main()
