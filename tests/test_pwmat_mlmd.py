import os
import numpy as np
import unittest
import dpdata


class TestSingleStep(unittest.TestCase):

    def setUp(self):
        self.LabeledSystem1 = dpdata.LabeledSystem(os.path.join('pwmat', 'OUT.MLMD'),\
        fmt='movement' )

    def test_mlmd(self) :

        self.assertEqual(self.LabeledSystem1['energies'], -0.2197270691E+03)
        self.assertEqual(self.LabeledSystem1.get_nframes(), 1)
        self.assertEqual(self.LabeledSystem1.get_natoms(), 5)
        self.assertEqual(self.LabeledSystem1.data['atom_names'], ['H', 'C'])
        self.assertEqual(self.LabeledSystem1.data['atom_numbs'], [4, 1])
    def test_cell(self) :
        fp = open('pwmat/mlmd_cell')
        cell = []
        for ii in fp :
            cell.append([float(jj) for jj in ii.split()])
        cell = np.array(cell)
        for ii in range(cell.shape[0]) :
            for jj in range(cell.shape[1]) :
                self.assertEqual(self.LabeledSystem1.data['cells'][0][ii][jj], cell[ii][jj])
        fp.close()
        
    def test_coord(self) :
        fp = open('pwmat/mlmd_coord')
        coord = []
        for ii in fp :
            coord.append([float(jj) for jj in ii.split()])
        coord = np.array(coord)
        for ii in range(coord.shape[0]) :
            for jj in range(coord.shape[1]) :
                self.assertEqual(self.LabeledSystem1.data['coords'][0][ii][jj], coord[ii][jj]*10.0)
        fp.close()
    def test_force(self) :
        fp = open('pwmat/mlmd_force')
        force = []
        for ii in fp :
            force.append([float(jj) for jj in ii.split()])
        force = np.array(force)
        for ii in range(force.shape[0]) :
            for jj in range(force.shape[1]) :
                self.assertEqual(self.LabeledSystem1.data['forces'][0][ii][jj], force[ii][jj])
        fp.close()



if __name__ == '__main__':
    unittest.main()
