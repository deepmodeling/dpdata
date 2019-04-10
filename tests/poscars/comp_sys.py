import numpy as np

class CompSys :

    def test_atom_numbs(self):
        self.assertEqual(self.system_1.data['atom_numbs'], 
                         self.system_2.data['atom_numbs'])

    def test_atom_names(self):
        self.assertEqual(self.system_1.data['atom_names'], 
                         self.system_2.data['atom_names'])

    def test_atom_types(self):
        self.assertEqual(self.system_1.data['atom_types'][0], 
                         self.system_2.data['atom_types'][0])
        self.assertEqual(self.system_1.data['atom_types'][1],
                         self.system_2.data['atom_types'][1])

    def test_orig(self):
        for d0 in range(3) :
            self.assertEqual(self.system_1.data['orig'][d0],
                             self.system_2.data['orig'][d0])

    def test_cell(self):
        for ii in range(3) :
            for jj in range(3) :
                self.assertAlmostEqual(self.system_1.data['cells'][0][ii][jj], 
                                       self.system_2.data['cells'][0][ii][jj], 
                                       places = 6,
                                       msg = 'cell[%d][%d] failed' % (ii,jj))

    def test_frame(self): 
        for ii in range(sum(self.system_1.data['atom_numbs'])) :
            for jj in range(3) :
                self.assertAlmostEqual(self.system_1.data['coords'][0][ii][jj], 
                                       self.system_2.data['coords'][0][ii][jj], 
                                       places = 6,
                                       msg = 'posis[%d][%d] failed' % (ii,jj))

