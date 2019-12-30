import os,shutil
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompLabeledSys, CompSys, IsPBC

class TestDeepmdLoadDumpComp(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp (self) :
        self.system_1 = dpdata.LabeledSystem('poscars/OUTCAR.h2o.md', 
                                             fmt = 'vasp/outcar')
        self.system_1.to_deepmd_npy('tmp.deepmd.npy', 
                                    prec = np.float64, 
                                    set_size = 2)        

        self.system_2 = dpdata.LabeledSystem('tmp.deepmd.npy', 
                                             fmt = 'deepmd/npy',
                                             type_map = ['O', 'H'])
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

    def tearDown(self) :
        if os.path.exists('tmp.deepmd.npy'):
            shutil.rmtree('tmp.deepmd.npy')


class TestDeepmdCompNoLabels(unittest.TestCase, CompSys, IsPBC) :
    def setUp (self) :
        self.system_1 = dpdata.System('poscars/POSCAR.h2o.md',
                                      fmt = 'vasp/poscar')
        self.system_1.to_deepmd_npy('tmp.deepmd.npy', 
                                    prec = np.float64, 
                                    set_size = 2)        
        self.system_2 = dpdata.System('tmp.deepmd.npy', 
                                      fmt = 'deepmd/npy',
                                      type_map = ['O', 'H'])
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

    def tearDown(self) :
        if os.path.exists('tmp.deepmd.npy'):
            shutil.rmtree('tmp.deepmd.npy')
        

class TestDeepmdCompNoLabels(unittest.TestCase, CompSys, IsPBC) :
    def setUp(self) :
        self.dir_name = 'tmp.deepmd.npy.nol'
        natoms = 3
        atom_names = ['O', 'H']
        atom_numbs = [1, 2]
        atom_types = np.array([0, 1, 1], dtype = np.int32)
        nframes = 11
        half_n = 6
        idx = [range(0, half_n), range(half_n, nframes)]
        os.makedirs(self.dir_name, exist_ok = True)
        os.makedirs(os.path.join(self.dir_name, 'set.000'), exist_ok = True)
        os.makedirs(os.path.join(self.dir_name, 'set.001'), exist_ok = True)
        np.savetxt(os.path.join(self.dir_name, 'type.raw'), atom_types, fmt = '%d')        
        
        coords = np.random.random([nframes, natoms, 3])
        cells = np.random.random([nframes, 3, 3])
        np.save(os.path.join(self.dir_name, 'set.000', 'coord.npy'), coords[idx[0]])
        np.save(os.path.join(self.dir_name, 'set.000', 'box.npy'),   cells [idx[0]])
        np.save(os.path.join(self.dir_name, 'set.001', 'coord.npy'), coords[idx[1]])
        np.save(os.path.join(self.dir_name, 'set.001', 'box.npy'),   cells [idx[1]])
        
        data = {
            'atom_names' : atom_names,
            'atom_types' : atom_types,
            'atom_numbs' : atom_numbs,
            'coords' : coords,
            'cells' : cells, 
            'orig' : np.zeros(3),
        }

        self.system_1 = dpdata.System(self.dir_name, fmt = 'deepmd/npy', type_map = ['O', 'H'])
        self.system_2 = dpdata.System()
        self.system_2.data = data

        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6


    def tearDown(self) :
        if os.path.exists(self.dir_name):
            shutil.rmtree(self.dir_name)


if __name__ == '__main__':
    unittest.main()
