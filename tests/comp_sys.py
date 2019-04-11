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

    def test_nframs(self):
        self.assertEqual(self.system_1.get_nframes(),
                         self.system_2.get_nframes())

    def test_cell(self):
        self.assertEqual(self.system_1.get_nframes(),
                         self.system_2.get_nframes())        
        for ff in range(self.system_1.get_nframes()) :
            for ii in range(3) :
                for jj in range(3) :
                    self.assertAlmostEqual(self.system_1.data['cells'][ff][ii][jj], 
                                           self.system_2.data['cells'][ff][ii][jj], 
                                           places = self.places,
                                           msg = 'cell[%d][%d][%d] failed' % (ff,ii,jj))

    def test_coord(self): 
        self.assertEqual(self.system_1.get_nframes(),
                         self.system_2.get_nframes())
        # think about direct coord
        tmp_cell = self.system_1.data['cells']
        tmp_cell = np.reshape(tmp_cell, [-1, 3])
        tmp_cell_norm = np.reshape(np.linalg.norm(tmp_cell, axis = 1), [-1, 3])
        for ff in range(self.system_1.get_nframes()) :
            for ii in range(sum(self.system_1.data['atom_numbs'])) :
                for jj in range(3) :
                    self.assertAlmostEqual(self.system_1.data['coords'][ff][ii][jj] / tmp_cell_norm[ff][jj], 
                                           self.system_2.data['coords'][ff][ii][jj] / tmp_cell_norm[ff][jj], 
                                           places = self.places,
                                           msg = 'coord[%d][%d][%d] failed' % (ff,ii,jj))


class CompLabeledSys (CompSys) :
    def test_energy(self) :
        self.assertEqual(self.system_1.get_nframes(),
                         self.system_2.get_nframes())
        for ff in range(self.system_1.get_nframes()) :
            self.assertAlmostEqual(self.system_1.data['energies'][ff], 
                                   self.system_2.data['energies'][ff], 
                                   places = self.places,
                                   msg = 'energies[%d] failed' % (ff))

    def test_force(self) :
        self.assertEqual(self.system_1.get_nframes(),
                         self.system_2.get_nframes())        
        for ff in range(self.system_1.get_nframes()) :
            for ii in range(self.system_1.data['forces'].shape[0]) :
                for jj in range(3) :
                    self.assertAlmostEqual(self.system_1.data['forces'][ff][ii][jj], 
                                           self.system_2.data['forces'][ff][ii][jj], 
                                           places = self.places,
                                           msg = 'forces[%d][%d][%d] failed' % (ff,ii,jj))
            
    def test_virial(self) :
        self.assertEqual(self.system_1.get_nframes(),
                         self.system_2.get_nframes())        
        for ff in range(self.system_1.get_nframes()) :
            for ii in range(3) :
                for jj in range(3) :
                    self.assertAlmostEqual(self.system_1.data['virials'][ff][ii][jj], 
                                           self.system_2.data['virials'][ff][ii][jj], 
                                           places = self.vir_places,
                                           msg = 'virials[%d][%d][%d] failed' % (ff,ii,jj))
            
                    
