import os
import numpy as np
import unittest
from context import dpdata
from dpdata.unit import LengthConversion

bohr2ang = LengthConversion("bohr", "angstrom").value()

class TestABACUSMD:

    def test_atom_names(self) :
        self.assertEqual(self.system_Si.data['atom_names'], ['Si'])
        #self.assertEqual(self.system_h2o.data['atom_names'], ['O','H'])

    def test_atom_numbs(self) :
        self.assertEqual(self.system_Si.data['atom_numbs'], [2])
        #self.assertEqual(self.system_h2o.data['atom_numbs'], [64,128])

    def test_atom_types(self) :
        ref_type = [0, 0]
        ref_type =  np.array(ref_type)
        for ii in range(ref_type.shape[0]) :
            self.assertEqual(self.system_Si.data['atom_types'][ii], ref_type[ii])

    def test_cell(self) :
        cell = bohr2ang * 10.2 * np.array([[0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
        for idx in range(np.shape(self.system_Si.data['cells'])[0]):
            for ii in range(cell.shape[0]) :
                for jj in range(cell.shape[1]) :
                    self.assertAlmostEqual(self.system_Si.data['cells'][idx][ii][jj], cell[ii][jj])

    def test_coord(self) :
        fp = open('abacus.md/Si_coord')
        coord = []
        for ii in fp :
            coord.append([float(jj) for jj in ii.split()])
        coord = np.array(coord)
        coord = coord.reshape([5, 2, 3])
        for ii in range(coord.shape[0]) :
            for jj in range(coord.shape[1]) :
                for kk in range(coord.shape[2]):
                    self.assertAlmostEqual(self.system_Si.data['coords'][ii][jj][kk], coord[ii][jj][kk])
        fp.close()

    def test_force(self) :
        fp = open('abacus.md/Si_force')
        force = []
        for ii in fp :
            force.append([float(jj) for jj in ii.split()])
        force = np.array(force)
        force = force.reshape([5, 2, 3])
        for ii in range(force.shape[0]) :
            for jj in range(force.shape[1]) :
                for kk in range(force.shape[2]):
                    self.assertAlmostEqual(self.system_Si.data['forces'][ii][jj][kk], force[ii][jj][kk])
        fp.close()

    def test_virial(self) :
        fp = open('abacus.md/Si_virial')
        virial = []
        for ii in fp :
            virial.append([float(jj) for jj in ii.split()])
        virial = np.array(virial)
        virial = virial.reshape([5, 3, 3])
        for ii in range(virial.shape[0]) :
            for jj in range(virial.shape[1]) :
                for kk in range(virial.shape[2]) :
                    self.assertAlmostEqual(self.system_Si.data['virials'][ii][jj][kk], virial[ii][jj][kk])
        fp.close()    

    def test_energy(self) :
        ref_energy = np.array([-211.77183266, -211.7739761 , -211.77713677, -211.78079673,
       -211.78511428])
        for ii in range(5):
            self.assertAlmostEqual(self.system_Si.data['energies'][ii], ref_energy[ii])


class TestABACUSMDLabeledOutput(unittest.TestCase, TestABACUSMD):

    def setUp(self):
        self.system_Si = dpdata.LabeledSystem('abacus.md',fmt='abacus/md')
        # self.system_h2o = dpdata.LabeledSystem('qe.scf/02.out',fmt='qe/pw/scf')


if __name__ == '__main__':
    unittest.main()