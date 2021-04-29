import os
import unittest
from context import dpdata
from rdkit import Chem

class TestBondOrderSystem(unittest.TestCase):

    def test_from_mol_file(self):
        syst = dpdata.BondOrderSystem("bond_order/CH3OH.mol", fmt='mol', type_map=['O','C','H'])
        self.assertEqual(syst.get_nframes(), 1)
        self.assertEqual(syst.get_nbonds(), 5)
        self.assertEqual(syst.get_natoms(), 6)
        self.assertEqual(syst['atom_names'], ['O','C','H'])
        self.assertAlmostEqual(syst['coords'][0][0][0], -0.3858)

    def test_regularize_formal_charges(self):
        non_regular = Chem.MolFromMolFile("bond_order/formal_charge.mol", removeHs=False)
        regular = dpdata.BondOrderSystem("bond_order/formal_charge.mol", fmt="mol")
        self.assertFalse(non_regular)
        self.assertTrue(isinstance(regular.rdkit_mol, Chem.rdchem.Mol))

    def test_read_other_format_without_bond_info(self):
        self.assertRaises(RuntimeError, dpdata.BondOrderSystem, "gromacs/1h.gro")
    
    def test_dump_to_deepmd_raw(self):
        syst = dpdata.BondOrderSystem("bond_order/formal_charge.mol", fmt="mol")
        syst.to_deepmd_raw("bond_order/deepmd_raw_test")
    
    def test_from_sdf_file(self):
        syst = dpdata.BondOrderSystem("bond_order/methane.sdf", type_map=['C','H'])
        self.assertEqual(syst.get_nframes(), 4)
        self.assertEqual(syst.get_nbonds(), 4)
        self.assertEqual(syst.get_natoms(), 5)
        self.assertEqual(syst['atom_names'], ['C','H'])
        self.assertAlmostEqual(syst['coords'][0][0][0], 0.0059)
        self.assertAlmostEqual(syst['coords'][1][0][0], 0.0043)
        self.assertAlmostEqual(syst['coords'][2][0][0], 0.0071)
        self.assertAlmostEqual(syst['coords'][3][0][0], 0.0032)
    
    def test_from_sdf_file_err(self):
        self.assertRaises(ValueError, dpdata.BondOrderSystem, "bond_order/methane_ethane.sdf")
