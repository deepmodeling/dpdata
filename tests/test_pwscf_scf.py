import os
import numpy as np
import unittest
from context import dpdata

class TestPWSCFSinglePointEnergy:
    def test_atom_names(self) :
        self.assertEqual(self.system.data['atom_names'], ['H','C'])
    def test_atom_numbs(self) :
        self.assertEqual(self.system.data['atom_numbs'], [4,1])
    def test_atom_types(self) :
        ref_type = [0,0,0,0,1]
        ref_type =  np.array(ref_type)
        for ii in range(ref_type.shape[0]) :
            self.assertEqual(self.system.data['atom_types'][ii], ref_type[ii])
    def test_cell(self) :
        cell = 10 * np.eye(3)
        for ii in range(cell.shape[0]) :
            for jj in range(cell.shape[1]) :
                self.assertEqual(self.system.data['cells'][0][ii][jj], cell[ii][jj])

    def test_coord(self) :
        fp = open('pwscf.scf/ch4_coord')
        coord = []
        for ii in fp :
            coord.append([float(jj) for jj in ii.split()])
        coord = np.array(coord)
        for ii in range(coord.shape[0]) :
            for jj in range(coord.shape[1]) :
                self.assertEqual(self.system.data['coords'][0][ii][jj], coord[ii][jj])
        fp.close()

    def test_force(self) :
        fp = open('pwscf.scf/ch4_force')
        force = []
        for ii in fp :
            force.append([float(jj) for jj in ii.split()])
        force = np.array(force)
        for ii in range(force.shape[0]) :
            for jj in range(force.shape[1]) :
                self.assertEqual(self.system.data['forces'][0][ii][jj], force[ii][jj])
        fp.close()

    def test_virial(self) :
        fp = open('pwscf.scf/ch4_virial')
        virial = []
        for ii in fp :
            virial.append([float(jj) for jj in ii.split()])
        virial = np.array(virial)
        for ii in range(virial.shape[0]) :
            for jj in range(virial.shape[1]) :
                self.assertEqual(self.system.data['virials'][0][ii][jj], virial[ii][jj])
        fp.close()

    def test_energy(self) :
        ref_energy = -219.74425946528794
        self.assertEqual(self.system.data['energies'][0], ref_energy)



class TestPWSCFLabeledOutput(unittest.TestCase, TestPWSCFSinglePointEnergy):

    def setUp(self):
        self.system = dpdata.LabeledSystem('pwscf.scf/01.out',fmt='pwscf/scf')

if __name__ == '__main__':
    unittest.main()

