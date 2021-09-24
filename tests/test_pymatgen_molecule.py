import os
import numpy as np
import unittest
from context import dpdata

class TestPOSCARCart(unittest.TestCase):
    
    def setUp(self): 
        self.system = dpdata.System()
        self.system.from_pymatgen_molecule(os.path.join('pymatgen', 'FA-001.xyz'))
        self.system.to("vasp/poscar", "FA-001.vasp")
        self.assertEqual(list(self.system["atom_types"]), [0, 1, 2, 1, 1, 2, 1, 1])

    def test_to_molecule(self):
        mols = self.system.to_pymatgen_molecule()
        self.assertEqual(len(mols), 1)

    def test_to_vasp(self):
        tmp_system = dpdata.System()
        tmp_system.from_vasp_poscar(os.path.join('pymatgen', 'mol2.vasp'))
        mols = tmp_system.to("pymatgen/molecule")
        struct = mols[-1].get_boxed_structure(10.0, 10.0, 10.0)
        struct.to(fmt = "poscar", filename = os.path.join('pymatgen', 'mol2-new.vasp'))


if __name__ == '__main__':
    unittest.main()
