import os
import numpy as np
import unittest, shutil
from context import dpdata
from dpdata.unit import LengthConversion

bohr2ang = LengthConversion("bohr", "angstrom").value()


class TestABACUSLabeledOutput(unittest.TestCase):
    def setUp(self):
        shutil.copy("abacus.scf/INPUT.ok", "abacus.scf/INPUT")
        self.system_ch4 = dpdata.LabeledSystem("abacus.scf", fmt="abacus/scf")
        # self.system_h2o = dpdata.LabeledSystem('qe.scf/02.out',fmt='qe/pw/scf')
        self.system_ch4_unlabeled = dpdata.System(
            "abacus.scf/STRU.ch4", fmt="abacus/stru"
        )

    def tearDown(self):
        if os.path.isfile("abacus.scf/INPUT"):
            os.remove("abacus.scf/INPUT")

    def test_atom_names(self):
        self.assertEqual(self.system_ch4.data["atom_names"], ["C", "H"])
        # self.assertEqual(self.system_h2o.data['atom_names'], ['O','H'])
        self.assertEqual(self.system_ch4_unlabeled.data["atom_names"], ["C", "H"])

    def test_atom_numbs(self):
        self.assertEqual(self.system_ch4.data["atom_numbs"], [1, 4])
        # self.assertEqual(self.system_h2o.data['atom_numbs'], [64,128])
        self.assertEqual(self.system_ch4_unlabeled.data["atom_numbs"], [1, 4])

    def test_atom_types(self):
        ref_type = [0, 1, 1, 1, 1]
        ref_type = np.array(ref_type)
        for ii in range(ref_type.shape[0]):
            self.assertEqual(self.system_ch4.data["atom_types"][ii], ref_type[ii])
            self.assertEqual(self.system_ch4_unlabeled["atom_types"][ii], ref_type[ii])
        # ref_type = [0]*64 + [1]*128
        # ref_type =  np.array(ref_type)
        # for ii in range(ref_type.shape[0]) :
        #     self.assertEqual(self.system_h2o.data['atom_types'][ii], ref_type[ii])

    def test_cell(self):
        # cell = 5.29177 * np.eye(3)
        cell = bohr2ang * 10 * np.eye(3)
        np.testing.assert_almost_equal(self.system_ch4.data["cells"][0], cell)
        np.testing.assert_almost_equal(self.system_ch4_unlabeled.data["cells"][0], cell)
        # fp = open('qe.scf/h2o_cell')
        # cell = []
        # for ii in fp :
        #     cell.append([float(jj) for jj in ii.split()])
        # cell = np.array(cell)
        # for ii in range(cell.shape[0]) :
        #     for jj in range(cell.shape[1]) :
        #         self.assertAlmostEqual(self.system_h2o.data['cells'][0][ii][jj], cell[ii][jj])
        # fp.close()

    def test_coord(self):
        with open("abacus.scf/ch4_coord") as fp:
            coord = []
            for ii in fp:
                coord.append([float(jj) for jj in ii.split()])
            coord = np.array(coord)
            np.testing.assert_almost_equal(
                self.system_ch4.data["coords"][0], coord, decimal=5
            )
            np.testing.assert_almost_equal(
                self.system_ch4_unlabeled.data["coords"][0], coord, decimal=5
            )

        # fp = open('qe.scf/h2o_coord')
        # coord = []
        # for ii in fp :
        #     coord.append([float(jj) for jj in ii.split()])
        # coord = np.array(coord)
        # for ii in range(coord.shape[0]) :
        #     for jj in range(coord.shape[1]) :
        #         self.assertAlmostEqual(self.system_h2o.data['coords'][0][ii][jj], coord[ii][jj])
        # fp.close()

    def test_force(self):
        with open("abacus.scf/ch4_force") as fp:
            force = []
            for ii in fp:
                force.append([float(jj) for jj in ii.split()])
            force = np.array(force)
            np.testing.assert_almost_equal(self.system_ch4.data["forces"][0], force)

        # fp = open('qe.scf/h2o_force')
        # force = []
        # for ii in fp :
        #     force.append([float(jj) for jj in ii.split()])
        # force = np.array(force)
        # for ii in range(force.shape[0]) :
        #     for jj in range(force.shape[1]) :
        #         self.assertAlmostEqual(self.system_h2o.data['forces'][0][ii][jj], force[ii][jj])
        # fp.close()

    def test_virial(self):
        with open("abacus.scf/ch4_virial") as fp:
            virial = []
            for ii in fp:
                virial.append([float(jj) for jj in ii.split()])
            virial = np.array(virial)
            np.testing.assert_almost_equal(
                self.system_ch4.data["virials"][0], virial, decimal=3
            )

        # fp = open('qe.scf/h2o_virial')
        # virial = []
        # for ii in fp :
        #     virial.append([float(jj) for jj in ii.split()])
        # virial = np.array(virial)
        # for ii in range(virial.shape[0]) :
        #     for jj in range(virial.shape[1]) :
        #         self.assertAlmostEqual(self.system_h2o.data['virials'][0][ii][jj], virial[ii][jj], places = 2)
        # fp.close()

    def test_energy(self):
        ref_energy = -219.64991404276591
        self.assertAlmostEqual(self.system_ch4.data["energies"][0], ref_energy)
        # ref_energy = -30007.651851226798
        # self.assertAlmostEqual(self.system_h2o.data['energies'][0], ref_energy)


class TestABACUSLabeledOutputFail(unittest.TestCase):
    def setUp(self):
        shutil.copy("abacus.scf/INPUT.fail", "abacus.scf/INPUT")
        self.system_ch4 = dpdata.LabeledSystem("abacus.scf", fmt="abacus/scf")
        # self.system_h2o = dpdata.LabeledSystem('qe.scf/02.out',fmt='qe/pw/scf')
        self.system_ch4_unlabeled = dpdata.System(
            "abacus.scf/STRU.ch4", fmt="abacus/stru"
        )

    def tearDown(self):
        if os.path.isfile("abacus.scf/INPUT"):
            os.remove("abacus.scf/INPUT")

    def test_return_zero(self):
        self.assertEqual(len(self.system_ch4), 0)


if __name__ == "__main__":
    unittest.main()
