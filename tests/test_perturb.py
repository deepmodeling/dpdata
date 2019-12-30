import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompSys, IsPBC

from unittest.mock import Mock
from unittest.mock import patch, MagicMock 

class NormalGenerator(object):
    def __init__(self):
        self.randn_generator = self.get_randn_generator()
        self.rand_generator = self.get_rand_generator()
    def randn(self,number):
        return next(self.randn_generator)
    def rand(self,number):
        return next(self.rand_generator)
    @staticmethod
    def get_randn_generator():
        data = np.asarray([
            [ 0.71878148, -2.20667426,  1.49373955],
            [-0.42728113,  1.43836059, -1.17553854],
            [-1.70793073, -0.39588759, -0.40880927],
            [ 0.17078291, -0.34856352,  1.04307936],
            [-0.99103413, -0.1886479,   0.13813131],
            [ 0.5839343,   1.04612646, -0.62631026],
            [ 0.9752889,   1.85932517, -0.47875828],
            [-0.23977172, -0.38373444, -0.04375488]])
        count = 0
        while True:
            yield data[count]
            count +=1
    
    @staticmethod        
    def get_rand_generator():
        yield np.asarray([0.23182233, 0.87106847, 0.68728511, 0.94180274, 0.92860453, 0.69191187])

class UniformGenerator(object):
    def __init__(self):
        self.randn_generator = self.get_randn_generator()
        self.rand_generator = self.get_rand_generator()
    def randn(self,number):
        return next(self.randn_generator)
    def rand(self,number):
        return next(self.rand_generator)

    @staticmethod
    def get_randn_generator():
        data = [[-0.19313281,  0.80194715,  0.14050915],
            [-1.47859926,  0.12921667, -0.17632456],
            [-0.60836805, -0.7700423,  -0.8386948 ],
            [-0.03236753,  0.36690245,  0.5041072 ],
            [-1.59366933,  0.37069227,  0.89608291],
            [ 0.18165617,  0.53875315, -0.42233955],
            [ 0.74052496,  1.26627555, -1.12094823],
            [-0.89610092, -1.44247021, -1.3502529 ]]
        yield np.asarray([0.0001,0.0001,0.0001]) # test for not using small vector 
        count = 0
        while True:
            yield data[count]
            count +=1

    @staticmethod        
    def get_rand_generator():
        data = np.asarray([[0.71263084], [0.61339295], 
            [0.22948181], [0.36087632],
            [0.17582222], [0.97926742],
            [0.84706761], [0.44495513]])
        
        yield np.asarray([0.34453551, 0.0618966,  0.9327273,  0.43013654, 0.88624993, 0.48827425])
        count =0
        while True:
            yield np.asarray(data[count])
            count+=1

class ConstGenerator(object):
    def __init__(self):
        self.randn_generator = self.get_randn_generator()
        self.rand_generator = self.get_rand_generator()
    def randn(self,number):
        return next(self.randn_generator)
    def rand(self,number):
        return next(self.rand_generator)

    @staticmethod
    def get_randn_generator():
        data = np.asarray([[ 0.95410606, -1.62338002, -2.05359934],
            [ 0.69213769, -1.26008667,  0.77970721],
            [-1.77926476, -0.39227219,  2.31677298],
            [ 0.08785233, -0.03966649, -0.45325656],
            [-0.53860887,  0.42536802, -0.46167309],
            [-0.26865791, -0.19901684, -2.51444768],
            [-0.31627314,  0.22076982, -0.36032225],
            [0.66731887, 1.2505806,  1.46112938]])
        yield np.asarray([0.0001,0.0001,0.0001]) # test for not using small vector 
        count = 0
        while True:
            yield data[count]
            count +=1

    @staticmethod        
    def get_rand_generator():        
        yield np.asarray([0.01525907, 0.68387374, 0.39768541, 0.55596047, 0.26557088, 0.60883073])

# %%
class TestPerturbNormal(unittest.TestCase, CompSys, IsPBC):
    @patch('numpy.random')
    def setUp (self, random_mock):
        random_mock.rand = NormalGenerator().rand
        random_mock.randn = NormalGenerator().randn
        system_1_origin = dpdata.System('poscars/POSCAR.SiC',fmt='vasp/poscar')
        self.system_1 = system_1_origin.perturb(1,0.05,0.6,'normal')
        self.system_2 = dpdata.System('poscars/POSCAR.SiC.normal',fmt='vasp/poscar')
        self.places = 6

class TestPerturbUniform(unittest.TestCase, CompSys, IsPBC):
    @patch('numpy.random')
    def setUp (self, random_mock) :
        random_mock.rand = UniformGenerator().rand
        random_mock.randn = UniformGenerator().randn
        system_1_origin = dpdata.System('poscars/POSCAR.SiC',fmt='vasp/poscar')
        self.system_1 = system_1_origin.perturb(1,0.05,0.6,'uniform')
        self.system_2 = dpdata.System('poscars/POSCAR.SiC.uniform',fmt='vasp/poscar')
        self.places = 6

class TestPerturbConst(unittest.TestCase, CompSys, IsPBC):
    @patch('numpy.random')
    def setUp (self, random_mock) :
        random_mock.rand = ConstGenerator().rand
        random_mock.randn = ConstGenerator().randn
        system_1_origin = dpdata.System('poscars/POSCAR.SiC',fmt='vasp/poscar')
        self.system_1 = system_1_origin.perturb(1,0.05,0.6,'const')
        self.system_2 = dpdata.System('poscars/POSCAR.SiC.const',fmt='vasp/poscar')
        self.places = 6

if __name__ == '__main__':
    unittest.main()