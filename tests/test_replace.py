import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompSys, IsPBC

from unittest.mock import Mock
from unittest.mock import patch, MagicMock 


class ConstGenerator(object):
    def __init__(self):
        self.choice_generator = self.get_choice_generator()
    def choice(self, a, size=None, replace=True, p=None):
        return next(self.choice_generator)

    @staticmethod
    def get_choice_generator():
        yield np.asarray([20,  6,  7, 22, 29,  2, 23, 10])

class TestReplace(unittest.TestCase, CompSys, IsPBC):
    @patch('numpy.random')
    def setUp(self, random_mock):
        random_mock.choice = ConstGenerator().choice
        self.system_1 = dpdata.System('poscars/POSCAR.P42nmc',fmt='vasp/poscar')
        self.system_1.replace('Hf', 'Zr', 8)
        # print(self.system_1.data)
        self.system_2 = dpdata.System('poscars/POSCAR.P42nmc.replace',fmt='vasp/poscar')
        # print(self.system_2.data)
        self.places = 6

# class TestReplicate123_not_change_origin(unittest.TestCase, CompSys, IsPBC):
#     def setUp (self) :
#        self.system_1 = dpdata.System('poscars/POSCAR.SiC',fmt='vasp/poscar')
#        self.system_1.replicate((1,2,3,))
#        self.system_2 = dpdata.System('poscars/POSCAR.SiC',fmt='vasp/poscar')
#        self.places = 6

if __name__ == '__main__':
    unittest.main()
