import os
import numpy as np
import unittest
from context import dpdata

class TestGromacsGro(unittest.TestCase):
    def test_read_file(self):
        system = dpdata.System('gromacs/1h.gro', type_map=['H', 'O'])
        self.assertTrue('H' in system['atom_names'])
        self.assertTrue('O' in system['atom_names'])
        self.assertEqual(system['atom_numbs'], [6, 3])
        for cc,ii in enumerate([1, 0, 0, 1, 0, 0, 1, 0, 0]):
            self.assertEqual(system['atom_types'][cc], ii)
        self.assertEqual(len(system['cells']), 1)
        self.assertEqual(len(system['coords']), 1)
        for ii in range(3):
            for jj in range(3):
                if ii != jj:
                    self.assertAlmostEqual(system['cells'][0][ii][jj], 0)
        self.assertAlmostEqual(system['cells'][0][0][0], 7.822838765564372)
        self.assertAlmostEqual(system['cells'][0][1][1], 7.353572647182051)
        self.assertAlmostEqual(system['cells'][0][2][2], 9.036518515423753)
        self.assertAlmostEqual(system['coords'][0][8][0], 7.43)
        self.assertAlmostEqual(system['coords'][0][8][1], 5.12)
        self.assertAlmostEqual(system['coords'][0][8][2], 3.36)

    def test_read_file_tri(self):
        system = dpdata.System('gromacs/1h.tri.gro', type_map=['H', 'O'])
        self.assertTrue('H' in system['atom_names'])
        self.assertTrue('O' in system['atom_names'])
        self.assertEqual(system['atom_numbs'], [6, 3])
        for cc,ii in enumerate([1, 0, 0, 1, 0, 0, 1, 0, 0]):
            self.assertEqual(system['atom_types'][cc], ii)
        self.assertEqual(len(system['cells']), 1)
        self.assertEqual(len(system['coords']), 1)
        count = 0
        for ii in range(3):
            for jj in range(3):
                if ii != jj:
                    self.assertAlmostEqual(system['cells'][0][ii][jj], count)
                    count += 1
        self.assertAlmostEqual(system['cells'][0][0][0], 7.822838765564372)
        self.assertAlmostEqual(system['cells'][0][1][1], 7.353572647182051)
        self.assertAlmostEqual(system['cells'][0][2][2], 9.036518515423753)
        self.assertAlmostEqual(system['coords'][0][8][0], 7.43)
        self.assertAlmostEqual(system['coords'][0][8][1], 5.12)
        self.assertAlmostEqual(system['coords'][0][8][2], 3.36)
        system.to('vasp/poscar', 'POSCAR')

class TestGromacsGroMultiFrames(unittest.TestCase):
    def test_read_file(self):
        system = dpdata.System('gromacs/multi_frames.gro', type_map=['H', 'O'])
        self.assertTrue('H' in system['atom_names'])
        self.assertTrue('O' in system['atom_names'])
        self.assertEqual(system['atom_numbs'], [6, 3])
        for cc,ii in enumerate([1, 0, 0, 1, 0, 0, 1, 0, 0]):
            self.assertEqual(system['atom_types'][cc], ii)
        self.assertEqual(len(system['cells']), 2)
        self.assertEqual(len(system['coords']), 2)
        for ii in range(3):
            for jj in range(3):
                if ii != jj:
                    self.assertAlmostEqual(system['cells'][0][ii][jj], 0) # frame no.1
                    self.assertAlmostEqual(system['cells'][1][ii][jj], 0) # frame no.2
        # frame no.1
        self.assertAlmostEqual(system['cells'][0][0][0], 7.822838765564372)
        self.assertAlmostEqual(system['cells'][0][1][1], 7.353572647182051)
        self.assertAlmostEqual(system['cells'][0][2][2], 9.036518515423753)
        self.assertAlmostEqual(system['coords'][0][8][0], 7.43)
        self.assertAlmostEqual(system['coords'][0][8][1], 5.12)
        self.assertAlmostEqual(system['coords'][0][8][2], 3.36)
        # frame no.2
        self.assertAlmostEqual(system['cells'][1][0][0], 7.822838765564372)
        self.assertAlmostEqual(system['cells'][1][1][1], 7.353572647182051)
        self.assertAlmostEqual(system['cells'][1][2][2], 9.036518515423753)
        self.assertAlmostEqual(system['coords'][1][8][0], 7.43)
        self.assertAlmostEqual(system['coords'][1][8][1], 5.12)
        self.assertAlmostEqual(system['coords'][1][8][2], 3.36)


class TestFormatAtomName(unittest.TestCase):
    def test_format_atom_name(self):
        system = dpdata.System("gromacs/case_for_format_atom_name.gro", fmt='gromacs/gro', type_map=['H','C','N','O','Cl'])
        self.assertEqual(system.formula, "H11C14N3O2Cl2")
    
    def test_no_format_atom_name(self):
        system = dpdata.System("gromacs/case_for_format_atom_name.gro", fmt='gromacs/gro', format_atom_name=False)
        atoms = ['CL1', 'H6', 'C4', 'C3', 'C6', 'C11', 'H10', 'C2', 'N3', 'C14',
                 'H7', 'H8', 'C13', 'H2', 'H1', 'H4', 'O2', 'H9', 'O1', 'N2', 'C9',
                 'H3', 'C5', 'H11', 'N1', 'C7', 'C10', 'CL2', 'H5', 'C1', 'C8','C12']
        for at in atoms:
            self.assertTrue(at in system['atom_names'])


class TestDumpGromacsGro(unittest.TestCase):
    def setUp(self):
        self.system = dpdata.System('gromacs/multi_frames.gro', type_map=['H', 'O'])
    
    def test_dump_single_frame(self):
        self.system.to_gromacs_gro('gromacs/tmp_1.gro', frame_idx=0)
        tmp = dpdata.System('gromacs/tmp_1.gro', type_map=['H', 'O'])
        self.assertEqual(tmp.get_nframes(), 1)    

    def test_dump_multi_frames(self):
        self.system.to_gromacs_gro('gromacs/tmp_2.gro')
        tmp = dpdata.System('gromacs/tmp_2.gro', type_map=['H', 'O'])
        self.assertEqual(tmp.get_nframes(), 2)
    
    def tearDown(self):
        if os.path.exists('gromacs/tmp_1.gro'):
            os.remove('gromacs/tmp_1.gro')
        if os.path.exists('gromacs/tmp_2.gro'):
            os.remove('gromacs/tmp_2.gro')
