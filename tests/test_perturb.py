import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompSys

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
            [-1.67097394,  1.00565478,  1.75499669],
            [ 0.40123822, -0.26832357,  0.51643976],
            [-0.70325165,  0.50154408,  1.98677553],
            [-0.13296745,  0.55434201,  1.54581177],
            [-0.44232103,  0.25854521, -0.43925698],
            [-0.86827367, -0.27676856,  0.44365476],
            [0.92187855, 1.09601508, 1.56535943],
            [-0.61866517, -0.96390108, -0.09593058]])
        count = 0
        while True:
            yield data[count]
            count +=1
    
    @staticmethod        
    def get_rand_generator():
        yield np.asarray([0.38114487,0.68166018,0.72880982,0.59149848,0.9529526,0.4403747 ])

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
        data = [[-0.6338632, -0.80482399, -0.51334222],
            [-1.18911393, -1.24521526, -1.34214803],
            [ 0.20640311, -0.52666381, -0.71929358],
            [ 0.51478689,  1.11296436, -0.61690281],
            [-0.77160527,  1.11681373,  0.31440837],
            [ 0.71732744, -1.56937797, -1.60701756],
            [ 0.6121493,  -0.06944649,  2.58409679],
            [-0.86858833, -1.44664272,  0.40813879]]
        yield np.asarray([0.0001,0.0001,0.0001]) # test for not using small vector 
        count = 0
        while True:
            yield data[count]
            count +=1

    @staticmethod        
    def get_rand_generator():
        data = np.asarray([[0.96216328], [0.50364801], 
            [0.10566744], [0.82063694],
            [0.00212345], [0.3384801],
            [0.75329772], [0.59874454]])
        
        yield np.asarray([0.24009826,0.9346383,0.54357125,0.05932848,0.96751515,0.78645846])
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
        data = [[-0.37656207, -0.50545915,  1.23769463],
            [-0.77863147,  0.05874498,  0.81611906],
            [0.557889,   0.36082455, 1.89695245],
            [ 1.07846313, -1.07323247,  0.18513304],
            [-0.19418326,  0.75757126, -1.90317556],
            [ 0.16438486,  1.40294282, -1.6183854 ],
            [-1.71033533, -0.73018772, -0.19946578],
            [-0.29332031, -2.00584426, -0.49622181]]
        yield np.asarray([0.0001,0.0001,0.0001]) # test for not using small vector 
        count = 0
        while True:
            yield data[count]
            count +=1

    @staticmethod        
    def get_rand_generator():        
        yield np.asarray([0.10497343,0.30020419,0.60592304, 0.85344648, 0.60080211, 0.22293347])

# %%
class TestPerturbNormal(unittest.TestCase, CompSys):
    @patch('numpy.random')
    def setUp (self, random_mock):
        random_mock.rand = NormalGenerator().rand
        random_mock.randn = NormalGenerator().randn
        system_1_origin = dpdata.System('poscars/POSCAR.SiC',fmt='vasp/poscar')
        self.system_1 = system_1_origin.perturb(1,0.05,0.1,'normal')
        self.system_2 = dpdata.System('poscars/POSCAR.SiC.normal',fmt='vasp/poscar')
        self.places = 6

class TestPerturbUniform(unittest.TestCase, CompSys):
    @patch('numpy.random')
    def setUp (self, random_mock) :
        random_mock.rand = UniformGenerator().rand
        random_mock.randn = UniformGenerator().randn
        system_1_origin = dpdata.System('poscars/POSCAR.SiC',fmt='vasp/poscar')
        self.system_1 = system_1_origin.perturb(1,0.05,0.1,'uniform')
        self.system_2 = dpdata.System('poscars/POSCAR.SiC.uniform',fmt='vasp/poscar')
        self.places = 6

class TestPerturbConst(unittest.TestCase, CompSys):
    @patch('numpy.random')
    def setUp (self, random_mock) :
        random_mock.rand = ConstGenerator().rand
        random_mock.randn = ConstGenerator().randn
        system_1_origin = dpdata.System('poscars/POSCAR.SiC',fmt='vasp/poscar')
        self.system_1 = system_1_origin.perturb(1,0.05,0.1,'const')
        self.system_2 = dpdata.System('poscars/POSCAR.SiC.const',fmt='vasp/poscar')
        self.places = 6

if __name__ == '__main__':
    unittest.main()