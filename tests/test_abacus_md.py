import os
import numpy as np
import unittest
from context import dpdata
from dpdata.unit import LengthConversion

bohr2ang = LengthConversion("bohr", "angstrom").value()

class TestABACUSMD:

    def test_atom_names(self) :
        self.assertEqual(self.system_water.data['atom_names'], ['H', 'O'])
        self.assertEqual(self.system_Si.data['atom_names'], ['Si'])

    def test_atom_numbs(self) :
        self.assertEqual(self.system_water.data['atom_numbs'], [2, 1])
        self.assertEqual(self.system_Si.data['atom_numbs'], [2])

    def test_atom_types(self) :
        ref_type = [0, 0, 1]
        ref_type =  np.array(ref_type)
        ref_type2 = np.array([0, 0])
        for ii in range(ref_type.shape[0]) :
            self.assertEqual(self.system_water.data['atom_types'][ii], ref_type[ii])
        for ii in range(ref_type2.shape[0]) :
            self.assertEqual(self.system_Si.data['atom_types'][ii], ref_type2[ii])

    def test_cell(self) :
        cell = bohr2ang * 28 * np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        cell2 = bohr2ang * 5.1 * np.array([[1, 1, 0], [1, 0, 1], [0, 1, 1]])
        for idx in range(np.shape(self.system_water.data['cells'])[0]):
            np.testing.assert_almost_equal(cell, self.system_water.data['cells'][idx], decimal = 5)
        for idx in range(np.shape(self.system_Si.data['cells'])[0]):
            np.testing.assert_almost_equal(self.system_Si.data['cells'][idx], cell2, decimal = 5)

    def test_coord(self) :
        with open('abacus.md/water_coord') as fp:
            coord = []
            for ii in fp :
                coord.append([float(jj) for jj in ii.split()])
            coord = np.array(coord)
            coord = coord.reshape([5, 3, 3])
            np.testing.assert_almost_equal(self.system_water.data['coords'], coord, decimal = 5)

        with open('abacus.md.nostress/Si_coord') as fp2:
            coord = []
            for ii in fp2 :
                coord.append([float(jj) for jj in ii.split()])
            coord = np.array(coord)
            coord = coord.reshape([4, 2, 3])
            np.testing.assert_almost_equal(self.system_Si.data['coords'], coord, decimal = 5)

    def test_force(self) :
        with open('abacus.md/water_force') as fp:
            force = []
            for ii in fp :
                force.append([float(jj) for jj in ii.split()])
            force = np.array(force)
            force = force.reshape([5, 3, 3])
            np.testing.assert_almost_equal(self.system_water.data['forces'], force, decimal=5)


        with open('abacus.md.nostress/Si_force') as fp2:
            force = []
            for ii in fp2 :
                force.append([float(jj) for jj in ii.split()])
            force = np.array(force)
            force = force.reshape([4, 2, 3])
            np.testing.assert_almost_equal(self.system_Si.data['forces'], force, decimal=5)


    def test_virial(self) :
        with open('abacus.md/water_virial') as fp:
            virial = []
            for ii in fp :
                virial.append([float(jj) for jj in ii.split()])
            virial = np.array(virial)
            virial = virial.reshape([5, 3, 3])
            np.testing.assert_almost_equal(self.system_water.data['virials'], virial, decimal=5) 

    def test_energy(self) :
        ref_energy = np.array([-466.69285117, -466.69929051, -466.69829826, -466.70364664,
       -466.6976083])
        ref_energy2 = np.array([-211.77184603, -211.78111966, -211.79681663, -211.79875524])
        np.testing.assert_almost_equal(self.system_water.data['energies'], ref_energy)
        np.testing.assert_almost_equal(self.system_Si.data['energies'], ref_energy2)

    def test_to_system(self):
        pp_file=["H.upf","O.upf"]
        numerical_orbital=["H.upf","O.upf"]
        numerical_descriptor="jle.orb"
        mass=[1.008,15.994]
        self.system_water.to(file_name="abacus.md/water_stru",fmt='abacus/stru',pp_file=pp_file,\
                numerical_orbital=numerical_orbital,numerical_descriptor=numerical_descriptor,\
                mass=mass)
        self.assertTrue(os.path.isfile('abacus.md/water_stru'))
        if os.path.isfile('abacus.md/water_stru'):
            with open('abacus.md/water_stru') as f:
                iline=0
                for iline,l in enumerate(f):
                    iline += 1
                self.assertEqual(iline,30)


class TestABACUSMDLabeledOutput(unittest.TestCase, TestABACUSMD):

    def setUp(self):
        self.system_water = dpdata.LabeledSystem('abacus.md',fmt='abacus/md') # system with stress
        self.system_Si = dpdata.LabeledSystem('abacus.md.nostress',fmt='abacus/md') # system without stress

    def tearDown(self):
        if os.path.isfile('abacus.md/water_stru'):
            os.remove('abacus.md/water_stru')

if __name__ == '__main__':
    unittest.main()
