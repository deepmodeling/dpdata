import os
import numpy as np
import unittest
from context import dpdata


class TestSIESTASinglePointEnergy:
    def test_atom_names(self):
        self.assertEqual(self.system.data['atom_names'], ['Si'])

    def test_atom_numbs(self):
        self.assertEqual(self.system.data['atom_numbs'], [64])

    def test_atom_types(self):
        ref_type = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        ref_type = np.array(ref_type)

        for ii in range(ref_type.shape[0]):
            # print(self.system.data['atom_types'][0][ii])
            self.assertAlmostEqual(self.system.data['atom_types'][ii], ref_type[ii])

    def test_cell(self):
        fp = open('siesta/aimd/cell')
        ref_cell = []
        for ii in fp:
            for jj in ii.split():
                ref_cell.append(float(jj))
        fp.close()

        cells = self.system.data['cells'].flatten()
        idx = 0
        for ii in range(len(cells)):
            self.assertAlmostEqual(cells[ii], float(ref_cell[ii]))

    def test_coord(self):
        fp = open('siesta/aimd/coord')
        ref_coord = []
        for ii in fp:
            for jj in ii.split():
                ref_coord.append(float(jj))
        fp.close()
        coords = self.system.data['coords'].flatten()
        for ii in range(len(coords)):
            self.assertAlmostEqual(coords[ii], float(ref_coord[ii]))

    def test_force(self):
        eV = 1
        angstrom = 1
        fp = open('siesta/aimd/force')
        ref_force = []
        for ii in fp:
            for jj in ii.split():
                ref_force.append(float(jj))

        fp.close()
        forces = self.system.data['forces'].flatten()
        for ii in range(len(forces)):
            self.assertAlmostEqual(forces[ii], float(ref_force[ii]))

    def test_viriale(self):
        toViri = 1
        vol = 1308.4268
        fp = open('siesta/aimd/virial')
        ref_virial = []
        for ii in fp:
            for jj in ii.split():
                ref_virial.append(float(jj))
        fp.close()
        virials = self.system.data['virials'].flatten()
        for ii in range(len(virials)):
            self.assertAlmostEqual(virials[ii], float(ref_virial[ii]) * toViri * vol)

    def test_energy(self):
        eV = 1
        fp = open('siesta/aimd/energy')
        ref_energy = []
        for ii in fp:
            for jj in ii.split():
                ref_energy.append(float(jj))
        fp.close()
        energy = self.system.data['energies']
        for ii in range(len(energy)):
            self.assertAlmostEqual(energy[ii], ref_energy[ii])


class TestAimdSIESTALabeledOutput(unittest.TestCase, TestSIESTASinglePointEnergy):

    def setUp(self):
        self.system = dpdata.LabeledSystem('siesta/aimd/output', fmt='siesta/aiMD_output')
        # self.system.data = dpdata.siesta.output.obtain_frame('siesta/siesta_output')


if __name__ == '__main__':
    unittest.main()
