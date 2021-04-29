import os
import unittest
from context import dpdata
from rdkit import Chem

class TestBondOrderSystem(unittest.TestCase):
    
    def setUp(self):
        pass

    def test_regularize_formal_charges(self):
        syst = dpdata.BondOrderSystem("bond_order/formal_charge.mol", fmt="mol")
        syst.to_mol_file("bond_order/tmp.mol")
        mol = Chem.MolFromMolFile("bond_order/tmp.mol")
        self.assertTrue(mol)

    def test_read_other_format(self):
        self.assertRaises(RuntimeError, dpdata.BondOrderSystem, "gromacs/1h.gro")
    
    def test_dump_to_deepmd_raw(self):
        syst = dpdata.BondOrderSystem("bond_order/tmp.mol", fmt="mol")
        syst.to_deepmd_raw("bond_order/deepmd_raw_test")
    
    def tearDown(self):
        os.remove("bond_order/tmp.mol")