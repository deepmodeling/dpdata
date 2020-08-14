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
                self.assertEqual(self.system.data['cells'][0][ii][jj], cell[ii][jj])


    def test_coord(self):
        coord = np.loadtxt('cp2k/ref_coord')
        for ii in range(coord.shape[0]) :
            for jj in range(coord.shape[1]) :
                self.assertEqual(self.system.data['coords'][0][ii][jj], coord[ii][jj])

    def test_force(self):
        #eV = 2.72113838565563E+01
        #angstrom = 5.29177208590000E-01
        force = np.loadtxt('cp2k/ref_force')
        for ii in range(force.shape[0]) :
            for jj in range(force.shape[1]) :
                self.assertEqual(self.system.data['forces'][0][ii][jj], force[ii][jj])

    def test_energy(self):
        #eV = 2.72113838565563E+01
        ref_energy = -48061.44424374075
        self.assertEqual(self.system.data['energies'][0], ref_energy)

    def test_virial(self):
        virial = np.loadtxt("cp2k/ref_virial")
        for ii in range(virial.shape[0]) :
            for jj in range(virial.shape[1]) :
                self.assertEqual(self.system.data['virials'][0][ii][jj], virial[ii][jj])




class TestCP2KLabeledOutput(unittest.TestCase, TestCP2KSinglePointEnergy):

    def setUp(self):
        self.system = dpdata.LabeledSystem('cp2k/cp2k_output', fmt = 'cp2k/output')

if __name__ == '__main__':
    unittest.main()

