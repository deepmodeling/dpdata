import os
import numpy as np
import unittest
from context import dpdata

class TestCP2KSinglePointEnergy:
    def test_atom_names(self):
        self.assertEqual(self.system.data['atom_names'], ['Fe','O'])
    def test_atom_numbs(self):
        self.assertEqual(self.system.data['atom_numbs'], [12,18])
    def test_atom_types(self):
        ref_type = np.loadtxt('cp2k/ref_type')
        for ii in range(ref_type.shape[0]) :
            self.assertEqual(self.system.data['atom_types'][ii], ref_type[ii])
    def test_cell(self):
        cell = np.loadtxt('cp2k/ref_cell')
        for ii in range(cell.shape[0]) :
            for jj in range(cell.shape[1]) :
                self.assertAlmostEqual(self.system.data['cells'][0][ii][jj], cell[ii][jj])


    def test_coord(self):
        coord = np.loadtxt('cp2k/ref_coord')
        for ii in range(coord.shape[0]) :
            for jj in range(coord.shape[1]) :
                self.assertAlmostEqual(self.system.data['coords'][0][ii][jj], coord[ii][jj])

    def test_force(self):
        #eV = 2.72113838565563E+01
        #angstrom = 5.29177208590000E-01
        force = np.loadtxt('cp2k/ref_force')
        for ii in range(force.shape[0]) :
            for jj in range(force.shape[1]) :
                self.assertAlmostEqual(self.system.data['forces'][0][ii][jj], force[ii][jj], places=6)

    def test_energy(self):
        #eV = 2.72113838565563E+01
        ref_energy = -48061.44846401638
        self.assertEqual(self.system.data['energies'][0], ref_energy)

    def test_virial(self):
        virial = np.loadtxt("cp2k/ref_virial")
        for ii in range(virial.shape[0]) :
            for jj in range(virial.shape[1]) :
                self.assertAlmostEqual(self.system.data['virials'][0][ii][jj], virial[ii][jj], places=6)




class TestCP2KLabeledOutput(unittest.TestCase, TestCP2KSinglePointEnergy):

    def setUp(self):
        self.system = dpdata.LabeledSystem('cp2k/cp2k_output', fmt = 'cp2k/output')

class TestNonCoveragedCP2KOutput:
    def setUp (self) :
        self.system = dpdata.LabeledSystem('cp2k/cp2k_nocon_output',
                                           fmt = 'cp2k/output')

    def test_atom_types(self) :
        self.assertEqual(self.system.data['atom_types'], [])

    def test_cells(self) :
        self.assertEqual(self.system.data['cells'], [])

    def test_coords(self) :
        self.assertEqual(self.system.data['coords'], [])

    def test_energies(self) :
        self.assertEqual(self.system.data['energies'], [])

    def test_forces(self) :
        self.assertEqual(self.system.data['forces'], [])

    def test_virials(self) :
        self.assertFalse('virials' in self.system.data)


if __name__ == '__main__':
    unittest.main()

