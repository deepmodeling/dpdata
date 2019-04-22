import os
import numpy as np
import unittest
from context import dpdata

class TestMSD (unittest.TestCase) :
    def setUp(self) :
        self.system = dpdata.System()
        self.system.data['atom_types'] = np.array([0,1])
        self.system.data['atom_names'] = ['O', 'H']
        nframes = 10
        cell_size = 5
        self.system.data['cells'] = np.tile(cell_size * np.eye(3),
                                            (nframes,1,1))
        self.system.data['coords'] = np.zeros([nframes, 2, 3])
        for ff in range(nframes) :
            self.system.data['coords'][ff][0] = 1.0 * ff * np.array([1,0,0])
            self.system.data['coords'][ff][1] = 2.0 * ff * np.array([1,0,0])
        self.system.data['coords'] = self.system.data['coords'] % cell_size
        
    def test_msd(self) :
        # print(self.system['atom_types'] == 0)
        msd = dpdata.md.msd.msd(self.system)
        msd0 = dpdata.md.msd.msd(self.system, self.system['atom_types'] == 0)
        msd1 = dpdata.md.msd.msd(self.system, self.system['atom_types'] == 1)
        # print(msd)
        ncomp = msd.shape[0]
        for ii in range(ncomp) :
            self.assertAlmostEqual(msd0[ii], 
                                   ii * ii,
                                   msg = 'msd0[%d]' % ii)
            self.assertAlmostEqual(msd1[ii], 
                                   ii * ii * 4,
                                   msg = 'msd1[%d]' % ii)
            self.assertAlmostEqual(msd[ii], 
                                   (msd0[ii]+msd1[ii]) * 0.5, 
                                   'msd[%d]' % ii)
