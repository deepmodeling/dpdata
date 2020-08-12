import numpy as np
import unittest
from context import dpdata


class TestFhi_aims_MD:
    def test_atom_names(self):
        self.assertEqual(self.system.data['atom_names'], ["C","H","O","N"])

    def test_atom_numbs(self):
        self.assertEqual(self.system.data['atom_numbs'], [32,36,8,4])

    def test_atom_types(self):
        ref_type = [0, 1, 1,]
        ref_type = np.array(ref_type)
        for ii in range(ref_type.shape[0]):
            self.assertAlmostEqual(self.system.data['atom_types'][ii], ref_type[ii])

    def test_cell(self):
        ref_cell=np.loadtxt('fhi_aims/ref_cell_md_m.txt')
        ref_cell=ref_cell.flatten()
        cells = self.system.data['cells'].flatten()
        idx = 0
        for ii in range(len(cells)):
            self.assertAlmostEqual(cells[ii], float(ref_cell[ii]))

    def test_coord(self):
        ref_coord=np.loadtxt('fhi_aims/ref_coord_md_m.txt')
        ref_coord=ref_coord.flatten()
        coords = self.system.data['coords'].flatten()
        for ii in range(len(coords)):
            self.assertAlmostEqual(coords[ii], float(ref_coord[ii]))

    def test_force(self):
        ref_force=np.loadtxt('fhi_aims/ref_force_md_m.txt')
        ref_force=ref_force.flatten()
        forces = self.system.data['forces'].flatten()
        for ii in range(len(forces)):
            self.assertAlmostEqual(forces[ii], float(ref_force[ii]))

    def test_energy(self):
        ref_energy=np.loadtxt('fhi_aims/ref_energy_md_m.txt')
        ref_energy=ref_energy.flatten()
        energy = self.system.data['energies']
        for ii in range(len(energy)):
            self.assertAlmostEqual(energy[ii], ref_energy[ii])


class TestFhi_aims_Output(unittest.TestCase, TestFhi_aims_MD):
    def setUp(self):
        self.system = dpdata.LabeledSystem('fhi_aims/output_multi_elements', fmt='fhi_aims/md')

if __name__ == '__main__':
    unittest.main()
