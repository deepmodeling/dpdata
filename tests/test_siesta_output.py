import os
import numpy as np
import unittest
from context import dpdata

class TestSIESTASinglePointEnergy:
    def test_atom_names(self) :
        self.assertEqual(self.system.data['atom_names'], ['H','C'])
    def test_atom_numbs(self) :
        self.assertEqual(self.system.data['atom_numbs'], [4, 1])
    def test_atom_types(self) :
        ref_type = [0,0,0,0,1]
        ref_type = np.array(ref_type)

        for ii in range(ref_type.shape[0]) :
            # print(self.system.data['atom_types'][0][ii])
            self.assertAlmostEqual(self.system.data['atom_types'][ii], ref_type[ii])

    def test_cell(self) :
        fp = open('siesta/scf/ref_cell')
        cell = []
        for ii in fp :
            for jj in ii.split():
                cell.append(float(jj))
        cell = np.array(cell)
        # print(cell)
        fp.close()
        res = self.system.data['cells'][0].flatten()
        for ii in range(len(cell)):
            self.assertAlmostEqual(res[ii], cell[ii])

    def test_coord(self) :
        fp = open('siesta/scf/ref_coord')
        coord = []
        for ii in fp:
            for jj in ii.split():
                coord.append(float(jj))
        coord = np.array(coord)
        fp.close()
        res = self.system.data['coords'][0].flatten()
        for ii in range(len(coord)) :
            self.assertAlmostEqual(res[ii], float(coord[ii]))

    def test_force(self) :
        eV = 1
        angstrom = 1
        fp = open('siesta/scf/ref_force')
        force = []
        for ii in fp:
            for jj in ii.split():
                force.append(float(jj))
        force = np.array(force)
        fp.close()
        res = self.system.data['forces'][0].flatten()
        for ii in range(len(force)):
            self.assertAlmostEqual(res[ii], float(force[ii]))

    def test_viriale(self) :
        toViri = 1
        fp = open('siesta/scf/ref_cell')
        cell = []
        for ii in fp:
            for jj in ii.split():
                cell.append(float(jj))
        cell = np.array(cell)
        cells = cell.reshape(3,3)
        fp.close()

        toVol = []
        for ii in cells:
            ### calucate vol
            toVol.append(np.linalg.det(cells))

        fp = open('siesta/scf/ref_virial')
        virial = []
        for ii in fp:
            for jj in ii.split():
                virial.append(float(jj) * toViri * toVol[0])
        virial = np.array(virial)
        fp.close()
        res = self.system.data['virials'][0].flatten()
        for ii in range(len(virial)):
            self.assertAlmostEqual(res[ii], float(virial[ii]))

    def test_energy(self) :
        eV = 1
        ref_energy = -219.1640
        self.assertAlmostEqual(self.system.data['energies'][0], ref_energy*eV)


class TestSIESTALabeledOutput(unittest.TestCase, TestSIESTASinglePointEnergy):

    def setUp(self):
        self.system = dpdata.LabeledSystem('siesta/scf/siesta_output', fmt = 'siesta/output')

if __name__ == '__main__':
    unittest.main()

