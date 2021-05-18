import os
import unittest
from context import dpdata
import glob
from rdkit import Chem
from rdkit.Chem import AllChem
import shutil
import numpy as np

class TestBondOrderSystem(unittest.TestCase):

    def test_from_rdkit_mol(self):
        mol = Chem.MolFromSmiles("CC")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMultipleConfs(mol, 10)
        system = dpdata.BondOrderSystem(rdkit_mol=mol)
        self.assertEqual(system.get_nframes(), 10)
        self.assertEqual(system.get_nbonds(), 7)

    def test_from_mol_file(self):
        syst = dpdata.BondOrderSystem("bond_order/CH3OH.mol", fmt='mol', type_map=['O','C','H'])
        self.assertEqual(syst.get_nframes(), 1)
        self.assertEqual(syst.get_nbonds(), 5)
        self.assertEqual(syst.get_natoms(), 6)
        self.assertEqual(syst['atom_names'], ['O','C','H'])
        self.assertAlmostEqual(syst['coords'][0][0][0], -0.3858)
    
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

    def test_regularize_formal_charges(self):
        non_regular = Chem.MolFromMolFile("bond_order/formal_charge.mol", removeHs=False)
        regular = dpdata.BondOrderSystem("bond_order/formal_charge.mol", fmt="mol")
        self.assertFalse(non_regular)
        self.assertTrue(isinstance(regular.rdkit_mol, Chem.rdchem.Mol))
    
    def test_formal_charge(self):
        names = ["C5H5-", "CH3CC-", "CH3NC", "CH3NH3+", "CH3NO2", "OCH3+",
                 "gly", "arg", "oxpy", "CH3OPO3_2-", "CH3PH3+", "CH3OAsO3_2-",
                 "CH3SH", "CH3_2SO", "CH3_2SO2", "CH3SO3-", "BOH4-"]
        charges = [-1, -1, 0, 1, 0, 1, 0, 1, 0, -2, 1, -2, 0, 0, 0, -1, -1]
        mols = [dpdata.BondOrderSystem(f"bond_order/{name}.mol") for name in names]
        self.assertEqual(charges, [mol.get_charge() for mol in mols])

    def test_read_other_format_without_bond_info(self):
        self.assertRaises(RuntimeError, dpdata.BondOrderSystem, "gromacs/1h.gro")
    
    def test_dump_to_deepmd_raw(self):
        syst = dpdata.BondOrderSystem("bond_order/methane.sdf", fmt="sdf")
        syst.to_deepmd_raw("bond_order/methane")
        formal_charges = list(np.loadtxt("bond_order/methane/formal_charges.raw"))
        self.assertTrue(formal_charges, [0 for _ in range(5)])
        bonds = np.loadtxt("bond_order/methane/bonds.raw")
        for bond_idx in range(4):
            for ii in range(3):
                self.assertEqual(syst['bonds'][bond_idx][ii], bonds[bond_idx][ii])
        shutil.rmtree("bond_order/methane")
    
    def test_dump_to_deepmd_npy(self):
        syst = dpdata.BondOrderSystem("bond_order/methane.sdf", fmt="sdf")
        syst.to_deepmd_npy("bond_order/methane")
        formal_charges = list(np.loadtxt("bond_order/methane/formal_charges.raw"))
        self.assertTrue(formal_charges, [0 for _ in range(5)])
        bonds = np.loadtxt("bond_order/methane/bonds.raw")
        for bond_idx in range(4):
            for ii in range(3):
                self.assertEqual(syst['bonds'][bond_idx][ii], bonds[bond_idx][ii])
        shutil.rmtree("bond_order/methane")
    
    def test_sanitize_mol_obabel(self):
        cnt = 0
        for sdf_file in glob.glob("bond_order/refined-set-ligands/obabel/*sdf"):
            syst = dpdata.BondOrderSystem(sdf_file, sanitize_level='high', verbose=False)
            if syst.rdkit_mol is None:
                cnt += 1
        self.assertEqual(cnt, 0)
    
    def test_sanitize_mol_origin(self):
        cnt = 0
        for sdf_file in glob.glob("bond_order/refined-set-ligands/origin/*sdf"):
            syst = dpdata.BondOrderSystem(sdf_file, sanitize_level='high', verbose=False)
            if syst.rdkit_mol is None:
                cnt += 1
        self.assertEqual(cnt, 0)
    
    def tearDown(self):
        if os.path.exists("tests/.cache"):
            shutil.rmtree("tests/.cache")
