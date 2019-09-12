import os
import numpy as np
import unittest
from context import dpdata

import numpy as np


class TestSetAtomTypes(unittest.TestCase):
    def setUp (self) :
        self.system_1 = dpdata.LabeledSystem('poscars/vasprun.h2o.md.10.xml')
        self.type_1 = self.system_1.get_atom_types()
        self.system_types = np.array([0,0,1,1,1,1]) 
        self.type_2 = self.system_1.map_atom_types(["H","C","O"])
        self.type_3 = self.system_1.map_atom_types({"H":2,"C":1,"O":3})

                
    def test_types_func_1(self):
        atom_types=np.array([2,2,0,0,0,0])
        atom_types_system_2=self.type_2
        atom_types_system_1=self.type_1
        for d0 in range(3) :
            self.assertEqual(atom_types[d0],
                             atom_types_system_2[d0])
        for d0 in range(3) :
            self.assertEqual(self.system_types[d0],
                             atom_types_system_1[d0])

    def test_types_func_2(self):
        atom_types=np.array([3,3,2,2,2,2])
        atom_types_system_3=self.type_3
        atom_types_system_1=self.type_1
        for d0 in range(3) :
            self.assertEqual(atom_types[d0],
                             atom_types_system_3[d0])
        for d0 in range(3) :
            self.assertEqual(self.system_types[d0],
                             atom_types_system_1[d0])

if __name__ == '__main__':
    unittest.main()
