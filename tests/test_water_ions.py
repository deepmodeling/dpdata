import os
import numpy as np
import unittest
from context import dpdata
try:
    import ase
    import ase.neighborlist    
    exist_ase=True
except:
    exist_ase=False

class TestIons(unittest.TestCase):
    
    def setUp(self): 
        self.system = dpdata.System()
        self.system.from_lammps_lmp(os.path.join('poscars', 'conf.waterion.lmp'), 
                                    type_map = ['O', 'H'])
        self.bonds = dpdata.md.water.compute_bonds(self.system.data['cells'][0],
                                                   self.system.data['coords'][0],
                                                   self.system.data['atom_types'])

    def test_ions_count(self) :
        no, noh, noh2, noh3, nh \
            = dpdata.md.water.find_ions(self.system.data['atom_types'], self.bonds)
        self.assertEqual(len(no), 0)
        self.assertEqual(len(noh), 1)
        self.assertEqual(len(noh2), 125)
        self.assertEqual(len(noh3), 1)
        self.assertEqual(len(nh), 0)
        self.assertEqual(noh[0], 0)
        self.assertEqual(noh3[0], 51)


@unittest.skipIf(not exist_ase, "skip TestAseComputeBond")
class TestAseComputeBond(unittest.TestCase):
    def setUp(self):
        self.system = dpdata.System()
        self.system.from_lammps_lmp(os.path.join('poscars', 'conf.waterion.lmp'), 
                                    type_map = ['O', 'H'])
        self.bonds = dpdata.md.water.compute_bonds_naive(self.system.data['cells'][0],
                                                         self.system.data['coords'][0],
                                                         self.system.data['atom_types'])
        self.bonds_ase = dpdata.md.water.compute_bonds_ase(self.system.data['cells'][0],
                                                           self.system.data['coords'][0],
                                                           self.system.data['atom_types'])

    def test_bond_identity(self):
        self.assertTrue(len(self.bonds), len(self.bonds_ase))
        for ii in range(len(self.bonds)):
            self.assertTrue(set(self.bonds[ii]) == set(self.bonds_ase[ii]))


        
if __name__ == '__main__':
    unittest.main()
    
