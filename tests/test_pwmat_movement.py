import os
import numpy as np
import unittest
import dpdata

class TestpwmatSinglePointEnergy:
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
        fp = open('pwmat/ref_cell')
        cell = []
        for ii in fp :
            cell.append([float(jj) for jj in ii.split()])
        cell = np.array(cell)
        for ii in range(cell.shape[0]) :
            for jj in range(cell.shape[1]) :
                self.assertEqual(self.system.data['cells'][0][ii][jj], cell[ii][jj])
        fp.close()


    def test_coord(self) :
        fp = open('pwmat/ref_coord')
        coord = []
        for ii in fp :
            coord.append([float(jj) for jj in ii.split()])
        coord = np.array(coord)
        for ii in range(coord.shape[0]) :
            for jj in range(coord.shape[1]) :
                self.assertEqual(self.system.data['coords'][0][ii][jj], coord[ii][jj]*10.0)
        fp.close()

    def test_force(self) :
        fp = open('pwmat/ref_force')
        force = []
        for ii in fp :
            force.append([float(jj) for jj in ii.split()])
        force = np.array(force)
        for ii in range(force.shape[0]) :
            for jj in range(force.shape[1]) :
                self.assertEqual(self.system.data['forces'][0][ii][jj], force[ii][jj])
        fp.close()

    def test_energy(self) :
        ref_energy = -0.2196929065E+03
        self.assertEqual(self.system.data['energies'][0], ref_energy)



class TestpwmatLabeledOutput(unittest.TestCase, TestpwmatSinglePointEnergy):

    def setUp(self):
        self.system = dpdata.LabeledSystem('pwmat/MOVEMENT', fmt = 'pwmat/MOVEMENT')

if __name__ == '__main__':
    unittest.main()

