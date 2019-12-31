import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompLabeledSys

class TestCp2kAimdOutput(unittest.TestCase, CompLabeledSys):
    def setUp (self) :
        self.system_1 = dpdata.LabeledSystem('cp2k/aimd',fmt='cp2k/aimd_output')
        self.system_2 = dpdata.LabeledSystem('cp2k/aimd/deepmd', fmt='deepmd/raw')
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4

if __name__ == '__main__':
    unittest.main()
