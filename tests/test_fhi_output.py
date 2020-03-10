import os
import numpy as np
import unittest
from context import dpdata

class TestFhi_aims:
    def test_atom_names(self) :
        self.assertEqual(self.system.data['atom_names'], ['B','N'])
    def test_atom_numbs(self) :
        self.assertEqual(self.system.data['atom_numbs'], [1, 1])
    def test_atom_types(self) :
        ref_type = [0,1]
        ref_type = np.array(ref_type)
        for ii in range(ref_type.shape[0]) :
            self.assertAlmostEqual(self.system.data['atom_types'][ii], ref_type[ii])

    def test_cell(self) :
        cell = np.loadtxt('fhi_aims/ref_cell.txt').flatten()
        res = self.system.data['cells'][0].flatten()
        for ii in range(len(cell)):
            self.assertAlmostEqual(res[ii], cell[ii])

    def test_coord(self) :
        coord = np.loadtxt('fhi_aims/ref_coord.txt').flatten()
        res = self.system.data['coords'][0].flatten()
        for ii in range(len(coord)) :
            self.assertAlmostEqual(res[ii], float(coord[ii]))

    def test_force(self) :
        force = np.loadtxt('fhi_aims/ref_force.txt').flatten()
        res = self.system.data['forces'][0].flatten()
        for ii in range(len(force)):
            self.assertAlmostEqual(res[ii], float(force[ii]))

   # def test_viriale(self) :
   #     toViri = 1
   #     fp = open('fhi_aims/ref_cell')
   #     cell = []
   #     for ii in fp:
   #         for jj in ii.split():
   #             cell.append(float(jj))
   #     cell = np.array(cell)
   #     cells = cell.reshape(3,3)
   #     fp.close()

   #     toVol = []
   #     for ii in cells:
   #         ### calucate vol
   #         toVol.append(np.linalg.det(cells))

   #     fp = open('fhi_aims/ref_virial')
   #     virial = []
   #     for ii in fp:
   #         for jj in ii.split():
   #             virial.append(float(jj) * toViri * toVol[0])
   #     virial = np.array(virial)
   #     fp.close()
   #     res = self.system.data['virials'][0].flatten()
   #     for ii in range(len(virial)):
   #         self.assertAlmostEqual(res[ii], float(virial[ii]))

    def test_energy(self) :
        ref_energy = -0.215215685892915E+04
        self.assertAlmostEqual(self.system.data['energies'][0], ref_energy,places = 6)


class TestFhiOutput(unittest.TestCase, TestFhi_aims):

    def setUp(self):
        self.system = dpdata.LabeledSystem('fhi_aims/out_scf', fmt = 'fhi_aims/scf')

if __name__ == '__main__':
    unittest.main()

