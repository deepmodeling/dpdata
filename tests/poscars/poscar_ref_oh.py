import numpy as np

class TestPOSCARoh :

    def test_atom_numbs(self):
        self.assertEqual(self.system.data['atom_numbs'], [1,1])

    def test_atom_names(self):
        self.assertEqual(self.system.data['atom_names'], ['O','H'])

    def test_atom_types(self):
        self.assertEqual(self.system.data['atom_types'][0], 0)
        self.assertEqual(self.system.data['atom_types'][1], 1)

    def test_orig(self):
        for d0 in range(3) :
            self.assertEqual(self.system.data['orig'][d0], 0)

    def test_cell(self):
        ovito_cell = np.array([[2.5243712, 0.0000000, 0.0000000], 
                               [1.2621856, 2.0430257, 0.0000000], 
                               [1.2874292, 0.7485898, 2.2254033]])
        for ii in range(3) :
            for jj in range(3) :
                self.assertAlmostEqual(self.system.data['cells'][0][ii][jj], 
                                       ovito_cell[ii][jj], 
                                       places = 6,
                                       msg = 'cell[%d][%d] failed' % (ii,jj))

    def test_frame(self): 
        ovito_posis = np.array([[0, 0, 0],
                                [1.2621856, 0.7018028, 0.5513885]])
        for ii in range(2) :
            for jj in range(3) :
                self.assertAlmostEqual(self.system.data['coords'][0][ii][jj], 
                                       ovito_posis[ii][jj], 
                                       places = 6,
                                       msg = 'posis[%d][%d] failed' % (ii,jj))
