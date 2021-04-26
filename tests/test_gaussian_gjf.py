import os
import unittest
from context import dpdata


class TestDumpGaussianGjf(unittest.TestCase):
    def setUp(self):
        self.system = dpdata.LabeledSystem('gaussian/methane.gaussianlog', 
                                           fmt = 'gaussian/log')
    
    def test_dump_to_gjf(self):
        self.system.to_gaussian_gjf("gaussian/tmp.gjf")
    
    def tearDown(self):
        if os.path.exists('gaussian/tmp.gjf'):
            os.remove('gaussian/tmp.gjf')