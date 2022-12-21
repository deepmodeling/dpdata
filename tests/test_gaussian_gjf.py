import unittest
import os

from context import dpdata


class TestGaussianGJF(unittest.TestCase):
    def setUp (self) :
        self.system = dpdata.LabeledSystem('poscars/OUTCAR.h2o.md', 
                                             fmt = 'vasp/outcar')
    
    def test_dump_gaussian_gjf(self):
        self.system.to_gaussian_gjf('tmp.gjf')
        os.remove('tmp.deepmd.hdf5')
