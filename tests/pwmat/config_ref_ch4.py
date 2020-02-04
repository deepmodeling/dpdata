import numpy as np

class Testconfigch4 :

    def test_atom_numbs(self):
        self.assertEqual(self.system.data['atom_numbs'], [4,1])

    def test_atom_names(self):
        self.assertEqual(self.system.data['atom_names'], ['H','C'])

    def test_atom_types(self):
        self.assertEqual(self.system.data['atom_types'][0], 0)
        self.assertEqual(self.system.data['atom_types'][1], 0)
        self.assertEqual(self.system.data['atom_types'][2], 0)
        self.assertEqual(self.system.data['atom_types'][3], 0)
        self.assertEqual(self.system.data['atom_types'][4], 1)

    def test_orig(self):
        for d0 in range(3) :
            self.assertEqual(self.system.data['orig'][d0], 0)

    def test_cell(self):
        ovito_cell = np.array([[10.000000, 0.0000000, 0.0000000], 
                               [0.0000000, 10.000000, 0.0000000], 
                               [0.0000000, 0.0000000, 10.000000]])
        for ii in range(3) :
            for jj in range(3) :
                self.assertAlmostEqual(self.system.data['cells'][0][ii][jj], 
                                       ovito_cell[ii][jj], 
                                       places = 6,
                                       msg = 'cell[%d][%d] failed' % (ii,jj))

    def test_frame(self): 
        ovito_posis = np.array([[0.53815434, 0.40686080, 0.36057301],
                                [0.39453966, 0.48032057, 0.43846884],
                                [0.55209243, 0.56545029, 0.44270874],
                                [0.52818530, 0.41641476, 0.53918266],
                                [0.50325059, 0.46725516, 0.44523234]])*10
        for ii in range(2) :
            for jj in range(3) :
                self.assertAlmostEqual(self.system.data['coords'][0][ii][jj], 
                                       ovito_posis[ii][jj], 
                                       places = 6,
                                       msg = 'posis[%d][%d] failed' % (ii,jj))
