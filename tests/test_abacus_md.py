import os
import numpy as np
import unittest
from context import dpdata
from dpdata.unit import LengthConversion

bohr2ang = LengthConversion("bohr", "angstrom").value()

class TestABACUSMD:

    def test_atom_names(self) :
        self.assertEqual(self.system_water.data['atom_names'], ['H', 'O'])

    def test_atom_numbs(self) :
        self.assertEqual(self.system_water.data['atom_numbs'], [2, 1])

    def test_atom_types(self) :
        ref_type = [0, 0, 1]
        ref_type =  np.array(ref_type)
        for ii in range(ref_type.shape[0]) :
            self.assertEqual(self.system_water.data['atom_types'][ii], ref_type[ii])

    def test_cell(self) :
        cell = bohr2ang * 28 * np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        for idx in range(np.shape(self.system_water.data['cells'])[0]):
            for ii in range(cell.shape[0]) :
                for jj in range(cell.shape[1]) :
                    self.assertAlmostEqual(self.system_water.data['cells'][idx][ii][jj], cell[ii][jj])

    def test_coord(self) :
        fp = open('abacus.md/water_coord')
        coord = []
        for ii in fp :
            coord.append([float(jj) for jj in ii.split()])
        coord = np.array(coord)
        coord = coord.reshape([5, 3, 3])
        for ii in range(coord.shape[0]) :
            for jj in range(coord.shape[1]) :
                for kk in range(coord.shape[2]):
                    self.assertAlmostEqual(self.system_water.data['coords'][ii][jj][kk], coord[ii][jj][kk])
        fp.close()

    def test_force(self) :
        fp = open('abacus.md/water_force')
        force = []
        for ii in fp :
            force.append([float(jj) for jj in ii.split()])
        force = np.array(force)
        force = force.reshape([5, 3, 3])
        for ii in range(force.shape[0]) :
            for jj in range(force.shape[1]) :
                for kk in range(force.shape[2]):
                    self.assertAlmostEqual(self.system_water.data['forces'][ii][jj][kk], force[ii][jj][kk])
        fp.close()

    def test_virial(self) :
        fp = open('abacus.md/water_virial')
        virial = []
        for ii in fp :
            virial.append([float(jj) for jj in ii.split()])
        virial = np.array(virial)
        virial = virial.reshape([5, 3, 3])
        for ii in range(virial.shape[0]) :
            for jj in range(virial.shape[1]) :
                for kk in range(virial.shape[2]) :
                    self.assertAlmostEqual(self.system_water.data['virials'][ii][jj][kk], virial[ii][jj][kk])
        fp.close()    

    def test_energy(self) :
        ref_energy = np.array([-466.69285117, -466.69929051, -466.69829826, -466.70364664,
       -466.6976083])
        for ii in range(5):
            self.assertAlmostEqual(self.system_water.data['energies'][ii], ref_energy[ii])


class TestABACUSMDLabeledOutput(unittest.TestCase, TestABACUSMD):

    def setUp(self):
        self.system_water = dpdata.LabeledSystem('abacus.md',fmt='abacus/md')


if __name__ == '__main__':
    unittest.main()