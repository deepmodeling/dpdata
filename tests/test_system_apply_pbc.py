import os
import numpy as np
import unittest
from context import dpdata

class TestPBC(unittest.TestCase) :
    def test_pbc(self) :
        nframes = 10
        natoms = 20
        data = {}
        data['coords'] = np.random.random([nframes, natoms, 3]) + [5, 5, 5]
        data['cells'] = np.tile(10 * np.eye(3), [nframes, 1, 1])
        data['cells'] += np.random.random([nframes, 3, 3])
        shift = 20 * (np.random.random([nframes, natoms, 3]) - 0.5)
        shift = shift.astype(int)
        bk_coord = np.copy(data['coords'])
        data['coords'] += np.matmul(shift, data['cells'])
        sys = dpdata.System()
        sys.data = data
        sys.apply_pbc()
        for ii in range(nframes) :
            for jj in range(natoms) :
                for dd in range(3) :
                    self.assertAlmostEqual(sys['coords'][ii][jj][dd],
                                           bk_coord[ii][jj][dd],
                                           msg = 'coord[%d][%d][%d] failed' % (ii,jj,dd))
                    
if __name__ == '__main__':
    unittest.main()

