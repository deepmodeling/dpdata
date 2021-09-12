import os
import numpy as np
import unittest
import dpdata
from dpdata.customsystem import CustomSystem

# from dpdata.molecule import Molecule
from pymatgen.core.structure import Molecule
from pymatgen.core.periodic_table import Element
from collections import Counter

class TestCustomSystem(unittest.TestCase):
    
    def setUp(self): 
        self.cs = CustomSystem('hybridcrystal/cubic-PbI3.vasp', fmt = 'vasp/poscar')
        self.mol = Molecule.from_file('hybridcrystal/FA-001.xyz')
        

    def test_addmol_from_obj(self):
        print("\n--------test_addmol_from_obj-------\n")
        self.cs.addmol(mol = self.mol, at = [0.5, 0.5, 0.5])
        self.cs.to('vasp/poscar', 'hybridcrystal/cubic-FAPbI3.vasp', frame_idx=-1)
        nmols = self.cs.get_nmols()
        self.assertEqual(nmols, 1)
        self.cs.addmol(mol = self.mol, at = [0.5, 0.5, 0.0])
        nmols = self.cs.get_nmols()
        self.assertEqual(nmols, 2)
        self.cs.to('vasp/poscar', 'hybridcrystal/cubic-FAPbI3-2mols.vasp', frame_idx=0)
        print(self.cs)

    def test_append(self):
        print("\n--------test_append-------\n")
        self.cs1 = CustomSystem("hybridcrystal/beta-PbI3.vasp", fmt = "vasp/poscar")
        self.cs1.addmol(mol = self.mol, at = [0.5, 0.5, 0.5])
        self.cs2 = CustomSystem("hybridcrystal/gamma-PbI3.vasp", fmt = "vasp/poscar")
        self.cs1.append(self.cs2)
        print(self.cs1)
        nframes = self.cs1.get_nframes()
        self.assertEqual(nframes, 2)



    def test_replicate(self):
        print("\n--------test_replicate-------\n")
        self.cs.addmol(mol = self.mol, at = [0.5, 0.5, 0.5])
        large_cs = self.cs.replicate((2,2,2))
        large_cs.to('vasp/poscar', 'hybridcrystal/cubic-FAPbI3-sc222.vasp', frame_idx=-1)
        print(large_cs)

    def test_rotate(self):
        print("\n--------test_rotate-------\n")
        self.cs.addmol(mol = self.mol, at = [0.5, 0.5, 0.5])
        large_cs = self.cs.replicate((2,2,2))
        large_cs.to('vasp/poscar', 'hybridcrystal/cubic-FAPbI3-sc222.vasp', frame_idx=-1)
        large_cs.rotate_mol_by_idx(0, axis = [0,1,0], angle = 90)
        large_cs.to('vasp/poscar', 'hybridcrystal/cubic-FAPbI3-sc222-rot0.vasp', frame_idx=-1)
        large_cs.rotate_mol()
        large_cs.to('vasp/poscar', 'hybridcrystal/cubic-FAPbI3-sc222-rot.vasp', frame_idx=-1)
        self.assertEqual(large_cs["atom_names"], ['Pb', 'I', 'C', 'N', 'H'])
        self.assertEqual(large_cs["atom_numbs"], [8, 24, 8, 16, 40])
        ref_type = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 4, 4, 4, 4, 4, 2, 3, 3, 4, 4, 4, 4, 4, 2, 3, 3, 4, 4, 4, 4, 4, 2, 3, 3, 4, 4, 4, 4, 4, 2, 3, 3, 4, 4, 4, 4, 4, 2, 3, 3, 4, 4, 4, 4, 4, 2, 3, 3, 4, 4, 4, 4, 4, 2, 3, 3, 4, 4, 4, 4, 4]
        self.assertEqual(list(large_cs["atom_types"]), ref_type)
        # print(large_cs)
        # print("atom_names = ")
        # print(large_cs["atom_names"])
        # print("atom_numbs = ")
        # print(large_cs["atom_numbs"])
        # print("atom_types = ")
        # print(large_cs["atom_types"])
    
    
    def test_addmol_dupElem(self):
        print("\n--------test_addmol_dupElem----------\n")
        self.cs.addmol('hybridcrystal/FA-001-I.xyz', at = [0.5, 0.5, 0.5])
        self.cs.to('vasp/poscar', 'hybridcrystal/cubic-CN2I5-PbI3.vasp', frame_idx=-1)
        print(self.cs)


if __name__ == '__main__':
    unittest.main()
